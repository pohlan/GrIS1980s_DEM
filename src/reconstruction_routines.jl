function prepare_problem(obs_file::String, imbie_mask::String, bedm_file::String, model_files::Vector{String})
    # load observations
    obs = ncread(obs_file, "Band1")

    # load masks
    I_no_ocean, I_obs = get_indices(obs, imbie_mask, bedm_file)

    # load model data
    Data_ice, nx, ny = read_model_data(;F,model_files,I_no_ocean)
    obs_flat_I = F.(obs[I_no_ocean][I_obs])
    return  Data_ice, obs_flat_I, I_no_ocean, I_obs, nx, ny
end

function solve_problem(Data_ice::Matrix{T}, obs_flat_I::Vector{T}, I_no_ocean::Vector{Int}, I_obs::Vector{Int}, nx::Int, ny::Int, r::Int, λ::Real, use_arpack::Bool) where T <: Real
    # centering model data
    Data_mean  = mean(Data_ice, dims=2)
    Data_ice  .= Data_ice .- Data_mean
    # centering observations with model mean
    x_data     = obs_flat_I .- Data_mean[I_obs]

    # compute SVD
    println("Computing the SVD..")
    r = min(r, size(Data_ice, 2)-1)
    if use_arpack
        B, _, _ = svds(Data_ice, nsv=r)
        U = B.U
        Σ = B.S
    else
        B = svd(Data_ice)
        @views U = B.U[:,1:r]
        @views Σ = B.S[1:r]
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
    dif[I_no_ocean[I_obs]] .= x_rec[I_obs] .- obs_flat_I
    err_mean         = mean(abs.(dif[I_no_ocean[I_obs]]))
    @printf("Mean absolute error: %1.1f m\n", err_mean)

    return x_rec, err_mean, dif
end

function solve_lsqfit(λ::Real, r::Int, gr::Int, imbie_mask::String, bedm_file::String, model_files::Vector{String}, obs_file::String, do_figures=false, use_arpack=false)
    Data_ice, obs_flat_I, I_no_ocean, I_obs, nx, ny = prepare_problem(obs_file, imbie_mask, bedm_file, model_files)
    x_rec, err_mean, dif                            = solve_problem(Data_ice, obs_flat_I, I_no_ocean, I_obs, nx, ny, r, λ, use_arpack)

    # retrieve matrix of reconstructed DEM
    dem_rec             = zeros(F, nx,ny)
    dem_rec[I_no_ocean] = x_rec
    # smoothing
    dem_rec = mapwindow(median, dem_rec, (5,5))
    dem_rec[dem_rec .== 0] .= no_data_value

    # save as nc file
    mkpath("output/")
    println("Saving file..")
    logλ = Int(round(log(10, λ)))
    filename = r < min(size(Data_ice)...)-100 ? "output/rec_lambda_1e$logλ"*"_g$gr"*"_r$r.nc" : "output/rec_lambda_1e$logλ"*"_g$gr"*"_fullsvd.nc"
    layername   = "surface"
    attributes  = Dict(layername => Dict{String, Any}("long_name" => "ice surface elevation",
                                                      "standard_name" => "surface_altitude",
                                                      "units" => "m")
                        )
    save_netcdf(filename, obs_file, [dem_rec], [layername], attributes)
    # plot and save difference between reconstruction and observations
    if do_figures
        save_netcdf("output/dif.nc", obs_file, [dif], ["dif"], Dict("dif" => Dict()))
        Plots.heatmap(reshape(dif,nx,ny)', cmap=:bwr, clims=(-200,200), cbar_title="[m]", title="reconstructed - observations", size=(700,900))
        Plots.savefig(filename[1:end-3]*".png")
    end

    return filename
end

function create_reconstructed_bedmachine(rec_file, bedmachine_file)
    # load datasets
    surfaceDEM = ncread(rec_file, "surface")
    bedDEM     = ncread(bedmachine_file, "bed")
    bedm_mask  = ncread(bedmachine_file, "mask")
    ice_mask   = (surfaceDEM .!= no_data_value) .&& (surfaceDEM .> bedDEM)

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
    new_mask = ones(eltype(bedm_mask), size(bedm_mask))
    new_mask[bedm_mask .== 0] .= 0 # ocean
    new_mask[ice_mask] .= 2
    new_mask[floating_mask]   .= 3

    # make sure DEM is zero everywhere on the ocean and equal to bed DEM on bedrock
    surfaceDEM[new_mask .== 0] .= no_data_value
    surfaceDEM[new_mask .== 1] .= bedDEM[new_mask .== 1]

    # calculate ice thickness
    h_ice = zeros(size(surfaceDEM)) .+ no_data_value
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

    return dest
end
