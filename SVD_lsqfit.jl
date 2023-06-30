using svd_IceSheetDEM, NetCDF, Glob, PyPlot, ImageFiltering, Statistics

const F = Float32 # Julia default is Float64 but that kills the process for the full training data set if r is too large

ARGS = ["--save",
        "--lambda", "1e5",
        "--obs", "data/aerodem_raw/aerodem_rm-filtered_geoid-corr_g1200.nc",
        "--r", "377",
        "--obs_band_name", "Band1"]
# 377 -> findfirst(cumsum(Σ)./sum(Σ).>0.9)

# retrieve command line arguments
parsed_args = parse_commandline(ARGS)
λ           = F(parsed_args["λ"])     # regularization
res         = parsed_args["res"]      # resolution
filepath    = parsed_args["train_folder"]
obs_file    = parsed_args["obs"]
band_name   = parsed_args["obs_band_name"]
r           = parsed_args["r"]

# load observations
obs_orig = ncread(obs_file, band_name)

# get indices for elevation < 400 m --> force aerodem there
ixx = 0 .< obs_orig .< 400
obs = copy(obs_orig)
obs[ixx] .= 0

# load masks
I_no_ocean, I_obs = get_indices(obs, res)

# load model data
model_files = glob(joinpath(filepath,"usurf_ex_gris_g" * res * "*_id_*YM.nc"))
Data_all, nx, ny = read_model_data(;F,model_files, which_files=1:3)

# solve least-squares problem
dem_rec, dif, err_mean = solve_lsqfit(Data_all, obs, λ, r, I_no_ocean, I_obs)

# set values at < 400 m elevation (ixx) equal to aerodem
dem_rec[ixx]   .= obs_orig[ixx]
# set very small values to zero as they are most likely ocean
dem_rec[dem_rec .< 5.] .= 0.

# smooth out (lots of jumps at data gaps where aerodem is enforced in adjacent pixels)
dem_smooth = mapwindow(median, dem_rec, (5,5))

if parsed_args["save"]
    # save as nc file
    mkpath("output/")
    println("Saving file..")
    logλ = Int(round(log(10, λ)))
    filename = "output/rec_lambda_1e$logλ"*"_g$res"*"_r$r.nc"
    dem_smooth = dem_smooth[:,end:-1:1]  # a bit of a hack; turn Greenland 'upside down' so that it is correct in the final file
    save_netcdf(dem_smooth; dest=filename, sample_path=obs_file)
    # plot and save difference between reconstruction and observations
    figure(figsize=(14,16))
    p = pcolormesh(reshape(dif,nx,ny)',cmap="bwr"); colorbar(label="[m]"); p.set_clim(-200,200)
    title("reconstructed - observations")
    savefig(filename[1:end-3]*".jpg")
end
