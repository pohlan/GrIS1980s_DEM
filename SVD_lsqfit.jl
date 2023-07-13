using svd_IceSheetDEM, NetCDF, Glob, PyPlot, ImageFiltering, Statistics

const F = Float32 # Julia default is Float64 but that kills the process for the full training data set if r is too large

function solve_it(ARGS, gr, imbie_mask, training_data_path, obs_file)
    # retrieve command line arguments
    parsed_args = parse_commandline(ARGS)
    λ           = F(parsed_args["λ"])     # regularization
    r           = parsed_args["r"]

    # load observations
    obs_orig = ncread(obs_file, "Band1")

    # get indices for elevation < 400 m --> force aerodem there
    ixx = 0 .< obs_orig .< 400
    obs = copy(obs_orig)
    obs[ixx] .= 0

    # load masks
    I_no_ocean, I_obs = get_indices(obs, imbie_mask)

    # load model data
    model_files = glob(training_data_path * "usurf_ex_gris_g$(gr)*YM.nc")
    Data_all, nx, ny = read_model_data(;F,model_files)

    # solve least-squares problem
    dem_rec, dif, err_mean = solve_lsqfit(Data_all, obs, λ, r, I_no_ocean, I_obs)

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
