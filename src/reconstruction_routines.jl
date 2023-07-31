using Printf, Statistics, LinearAlgebra, TSVD, HDF5, ImageFiltering, PyPlot, NetCDF
import ArchGDAL as AG

function prepare_problem(obs_file, imbie_mask, model_files, F)
    # load observations
    obs_orig = ncread(obs_file, "Band1")

    # get indices for elevation < 400 m --> force aerodem there
    ixx = 0 .< obs_orig .< 400
    obs = copy(obs_orig)
    obs[ixx] .= 0

    # load masks
    I_no_ocean, I_obs = get_indices(obs, imbie_mask)

    # load model data
    Data_all, nx, ny = read_model_data(;F,model_files)

    return  Data_all, obs, I_no_ocean, I_obs, nx, ny, ixx, obs_orig
end

function solve_problem(Data_all, obs, I_no_ocean, I_obs, nx, ny, r, λ, F)
    # centering model data
    Data       = Data_all[I_no_ocean, :]  # remove cells where there is ocean, saves half of the space
    Data_mean  = mean(Data, dims=2)
    Data_centr = Data .- Data_mean
    # centering observations with model mean
    obs_flat   = obs[I_no_ocean]
    x_data     = obs_flat .- Data_mean

    # compute SVD
    println("Computing the SVD..")
    if r < min(size(Data_centr)...)-100  # the tsvd algorithm doesn't give good results for a full or close to full svd (https://github.com/JuliaLinearAlgebra/TSVD.jl/issues/28)
        U, Σ, _ = tsvd(Data_centr, r)
    else
        U, Σ, _ = svd(Data_centr)
    end

    # solve the lsqfit problem
    println("Solving the least squares problem..")
    UΣ            = U*diagm(Σ)
    U_A, Σ_A, V_A = svd(UΣ[I_obs,:])
    D             = transpose(diagm(Σ_A))*diagm(Σ_A) + λ*I
    v_rec         = V_A * D^(-1) * transpose(diagm(Σ_A)) * transpose(U_A) * x_data[I_obs]
    x_rec         = U*diagm(Σ)*v_rec .+ Data_mean

    # calculate error and print
    dif                     = zeros(F, nx,ny)
    dif[I_no_ocean[I_obs]] .= x_rec[I_obs] .- obs_flat[I_obs]
    err_mean         = mean(abs.(dif[I_no_ocean[I_obs]]))
    @printf("Mean absolute error: %1.1f m\n", err_mean)

    return x_rec, err_mean, dif
end

function solve_lsqfit(F, λ, r, gr, imbie_mask, model_files, obs_file)
    Data_all, obs, I_no_ocean, I_obs, nx, ny, ixx, obs_orig = prepare_problem(obs_file, imbie_mask, model_files, F)
    x_rec, err_mean, dif                                    = solve_problem(Data_all, obs, I_no_ocean, I_obs, nx, ny, r, λ, F)

    # retrieve matrix of reconstructed DEM
    dem_rec             = zeros(F, nx,ny)
    dem_rec[I_no_ocean] = x_rec

    # set values at < 400 m elevation (ixx) equal to aerodem
    dem_rec[ixx]   .= obs_orig[ixx]
    # set very small values to zero as they are most likely ocean
    dem_rec[dem_rec .< 5.] .= 0.

    # smooth out (lots of jumps at data gaps where aerodem is enforced in adjacent pixels)
    dem_smooth = mapwindow(median, dem_rec, (5,5))

    # save as nc file
    mkpath("output/")
    println("Saving file..")
    logλ = Int(round(log(10, λ)))
    filename = "output/rec_lambda_1e$logλ"*"_g$gr"*"_r$r.nc"
    dem_smooth = dem_smooth[:,end:-1:1]  # a bit of a hack; turn Greenland 'upside down' so that it is correct in the final file
    save_netcdf(dem_smooth; dest=filename, sample_path=obs_file)

    # plot and save difference between reconstruction and observations
    figure(figsize=(14,16))
    p = pcolormesh(reshape(dif,nx,ny)',cmap="bwr"); colorbar(label="[m]"); p.set_clim(-200,200)
    title("reconstructed - observations")
    savefig(filename[1:end-3]*".jpg")

    return filename
end

function create_reconstructed_bedmachine(obs_file, rec_file, bedmachine_file, template_path)
    # load datasets
    aeroDEM    = shortread(obs_file)
    surfaceDEM = shortread(rec_file)
    fid = h5open(bedmachine_file)
    bedDEM = read(fid["bed"])[:,end:-1:1]
    bedm_mask = read(fid["mask"])[:,end:-1:1]
    close(fid)
    ice_mask   = (surfaceDEM .> 0.0) .&& (surfaceDEM .> bedDEM)

    # calculate floating mask
    ρw            = 1030
    ρi            = 917
    Pw            = - ρw * bedDEM
    Pi            = ρi * (surfaceDEM - bedDEM)
    floating_mask = (Pw .> Pi) .&& (ice_mask)  # floating where water pressure > ice pressure at the bed

    # calculate mask
    new_mask = zeros(Int, size(bedm_mask))
    new_mask[ice_mask] .= 2
    new_mask[bedm_mask .== 1] .= 1 # bedrock
    new_mask[floating_mask]   .= 3

    # calculate ice thickness
    h_ice = zeros(size(surfaceDEM))
    h_ice[ice_mask] .= surfaceDEM[ice_mask] - bedDEM[ice_mask]
    h_ice[floating_mask] .= surfaceDEM[floating_mask] ./  (1-ρi/ρw)

    # save to netcdf file
    dest = "output/bedmachine1980_reconstructed.nc"
    template_dataset = AG.read(template_path)
        AG.create(
            dest,
            driver = AG.getdriver(template_dataset),
            width  = AG.width(template_dataset),
            height = AG.height(template_dataset),
            nbands = 4,
            dtype  = Float32,
            options = ["-co", "COMPRESS=DEFLATE", "-co", "ZLEVEL=6"] # reduces the file size
        ) do raster
            AG.write!(raster, surfaceDEM, 1)
            AG.write!(raster, bedDEM, 2)
            AG.write!(raster, h_ice, 3)
            AG.write!(raster, new_mask, 4)
            AG.setgeotransform!(raster, AG.getgeotransform(template_dataset))
            AG.setproj!(raster, AG.getproj(template_dataset))
        end

    # couldn't figure out how to name a layer with ArchGDAL
    for (n, varb) in enumerate(["surface", "bed", "thickness", "mask"])
        run(`ncrename -v $("Band"*string(n)),$varb $dest`)
    end
    return
end
