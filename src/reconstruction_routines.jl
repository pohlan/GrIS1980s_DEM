using Printf, Statistics, LinearAlgebra, TSVD, ImageFiltering, Plots, NetCDF, NCDatasets
import ArchGDAL as AG

function prepare_problem(obs_file, imbie_mask, model_files, F)
    # load observations
    obs = ncread(obs_file, "Band1")

    # get indices for elevation < 400 m --> force aerodem there
    ixx = 0 .< obs .< 400
    obs_force = obs[ixx]
    obs[ixx] .= 0

    # load masks
    I_no_ocean, I_obs = get_indices(obs, imbie_mask)

    # load model data
    Data_ice, nx, ny = read_model_data(;F,model_files,I_no_ocean)
    obs_flat_I = obs[I_no_ocean][I_obs]
    return  Data_ice, obs_flat_I, I_no_ocean, I_obs, nx, ny, ixx, obs_force
end

function solve_problem(Data_ice, obs_flat_I, I_no_ocean, I_obs, nx, ny, r, λ, F)
    # centering model data
    Data_mean   = mean(Data_ice, dims=2)
    Data_centr  = Data_ice .- Data_mean
    # centering observations with model mean
    x_data     = obs_flat_I .- Data_mean[I_obs]

    # compute SVD
    println("Computing the SVD..")
    B = svd(Data_centr)
    if r < min(size(Data_centr)...)-100  # the tsvd algorithm doesn't give good results for a full or close to full svd (https://github.com/JuliaLinearAlgebra/TSVD.jl/issues/28)
        # U, Σ, _ = tsvd(Data_centr, r)
        @views U = B.U[:,1:r]
        @views Σ = B.S[1:r]
    else
        U = B.U
        Σ = B.S
    end

    # solve the lsqfit problem
    println("Solving the least squares problem..")
    UΣ            = U*diagm(Σ)
    U_A, Σ_A, V_A = svd(UΣ[I_obs,:])
    D             = transpose(diagm(Σ_A))*diagm(Σ_A) + λ*I
    v_rec         = V_A * D^(-1) * transpose(diagm(Σ_A)) * transpose(U_A) * x_data
    x_rec         = U*diagm(Σ)*v_rec .+ Data_mean

    # calculate error and print
    dif                     = zeros(F, nx,ny)
    dif[I_no_ocean[I_obs]] .= obs_flat_I .- x_rec[I_obs]
    err_mean         = mean(abs.(dif[I_no_ocean[I_obs]]))
    @printf("Mean absolute error: %1.1f m\n", err_mean)

    return x_rec, err_mean, dif
end

function solve_lsqfit(F, λ, r, gr, imbie_mask, model_files, obs_file, do_figure=false)
    Data_ice, obs_flat_I, I_no_ocean, I_obs, nx, ny, ixx, obs_force = prepare_problem(obs_file, imbie_mask, model_files, F)
    x_rec, err_mean, dif                                    = solve_problem(Data_ice, obs_flat_I, I_no_ocean, I_obs, nx, ny, r, λ, F)

    # retrieve matrix of reconstructed DEM
    dem_rec             = zeros(F, nx,ny)
    dem_rec[I_no_ocean] = x_rec

    # set values at < 400 m elevation (ixx) equal to aerodem
    dem_rec[ixx] = obs_force
    # set very small values to zero as they are most likely ocean
    dem_rec[dem_rec .< 5.] .= 0.

    # smooth out (lots of jumps at data gaps where aerodem is enforced in adjacent pixels)
    dem_smooth = mapwindow(median, dem_rec, (5,5))

    # save as nc file
    mkpath("output/")
    println("Saving file..")
    logλ = Int(round(log(10, λ)))
    filename = r < min(size(Data_ice)...)-100 ? "output/rec_lambda_1e$logλ"*"_g$gr"*"_r$r.nc" : "output/rec_lambda_1e$logλ"*"_g$gr"*"_fullsvd.nc"
    layername   = "surface"
    attributes  = Dict(layername => Dict("long_name" => "ice surface elevation",
                                         "standard_name" => "surface_altitude",
                                         "units" => "m")
                        )
    save_netcdf(filename, obs_file, [dem_smooth], [layername], attributes)
    # plot and save difference between reconstruction and observations
    if do_figure
        heatmap(reshape(dif,nx,ny)', cmap=:bwr, clims=(-200,200), cbar_title="[m]", title="reconstructed - observations", size=(700,900))
        savefig(filename[1:end-3]*".pdf")
    end

    return filename
end

function create_reconstructed_bedmachine(rec_file, bedmachine_file)
    # load datasets
    surfaceDEM = ncread(rec_file, "surface")
    bedDEM     = ncread(bedmachine_file, "bed")
    bedm_mask  = ncread(bedmachine_file, "mask")
    ice_mask   = (surfaceDEM .> 0.0) .&& (surfaceDEM .> bedDEM)

    # retrieve grid size
    x          = ncread(rec_file, "x")
    gr         = Int(x[2] - x[1])

    # calculate floating mask
    ρw            = 1030
    ρi            = 917
    Pw            = - ρw * bedDEM
    Pi            = ρi * (surfaceDEM - bedDEM)
    floating_mask = (Pw .> Pi) .&& (ice_mask)  # floating where water pressure > ice pressure at the bed

    # calculate mask
    new_mask = zeros(typeof(bedm_mask[1,1]), size(bedm_mask))
    new_mask[ice_mask] .= 2
    new_mask[bedm_mask .== 1] .= 1 # bedrock
    new_mask[floating_mask]   .= 3

    # calculate ice thickness
    h_ice = zeros(size(surfaceDEM))
    h_ice[ice_mask] .= surfaceDEM[ice_mask] - bedDEM[ice_mask]
    h_ice[floating_mask] .= surfaceDEM[floating_mask] ./  (1-ρi/ρw)

    # save to netcdf file
    dest        = "output/bedmachine1980_reconstructed_g$(gr).nc"
    layers      = [surfaceDEM, bedDEM, h_ice, new_mask]
    layernames  = ["surface", "bed", "thickness", "mask"]
    template    = NCDataset(bedmachine_file)
    attributes  = get_attr(template, layernames)
    # overwrite some attributes
    sources_rec = Dict("surface"   => "svd reconstruction",
                       "bed"       => "Bedmachine-v5: Morlighem et al. (2022). IceBridge BedMachine Greenland, Version 5. Boulder, Colorado USA. NASA National Snow and Ice Data Center Distributed Active Archive Center. https://doi.org/10.5067/GMEVBWFLWA7X; projected on new grid with gdalwarp",
                       "thickness" => "computed from surface and bed",
                       "mask"      => "bedrock from Morlighem et al. (2022); ice, floating and ocean computed from surface and bed elevation"
                       )
    for l in layernames
        attributes[l]["source"] = sources_rec[l]
    end
    attributes["mask"]["long_name"] = "mask (0 = ocean, 1 = ice-free land, 2 = grounded ice, 3 = floating ice)"
    save_netcdf(dest, bedmachine_file, layers, layernames, attributes)

    return
end
