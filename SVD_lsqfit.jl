using svd_IceSheetDEM, NetCDF, Statistics, LinearAlgebra, Glob, PyPlot, Printf, TSVD, ImageFiltering

ARGS = ["--save", "--lambda", "1e5", "--obs", "data/aerodem_raw/aerodem_rm-filtered_geoid-corr_g1200.nc", "--r", "377",
        "--obs_band_name", "Band1"]
# 377 -> findfirst(cumsum(Σ)./sum(Σ).>0.9)

# retrieve command line arguments
parsed_args = parse_commandline(ARGS)
λ           = parsed_args["λ"]        # regularization
res         = parsed_args["res"]      # resolution
filepath    = parsed_args["train_folder"]
obs_file    = parsed_args["obs"]
band_name   = parsed_args["obs_band_name"]
r           = parsed_args["r"]

F        = Float32 # Julia default is Float64 but that kills the process for the full training data set if r is too large

# load observations
obs_orig = ncread(obs_file, band_name)

# get indices for elevation < 400 m --> force aerodem there
ixx = 0 .< obs_orig .< 400
obs = copy(obs_orig)
obs[ixx] .= 0

# load masks
I_no_ocean, I_data, I_intr = get_indices(obs, res)

# load model data
model_files = glob(joinpath(filepath,"usurf_ex_gris_g" * res * "*_id_*YM.nc"))
Data_all, nx, ny = read_model_data(;F,model_files)

# centering model data
Data       = Data_all[I_no_ocean, :]   # remove cells where there is ocean, saves half of the space
Data_mean  = mean(Data, dims=2)
Data_centr = Data .- Data_mean

# centering observations with model mean
obs_flat   = F.(reshape(obs, nx * ny, 1)[I_no_ocean])
x_data     = obs_flat .- Data_mean

# compute SVD
println("Computing the SVD..")
U, Σ, V = svd(Data_centr)

# solve the lsqfit problem
println("Solving the least squares problem..")
UΣ            = U*diagm(Σ)
U_A, Σ_A, V_A = svd(UΣ[I_data,:])
D             = transpose(diagm(Σ_A))*diagm(Σ_A) + λ*I
v_rec         = V_A * D^(-1) * transpose(diagm(Σ_A)) * transpose(U_A) * x_data[I_data]
x_rec         = U*diagm(Σ)*v_rec .+ Data_mean

# calculate error and print
dif                      = zeros(nx*ny)
dif[I_no_ocean[I_data]] .= (x_rec .- obs_flat)[I_data]
if !startswith(obs_file, filepath*"aerodem")
    dif[I_no_ocean[I_intr]] .= (x_rec .- obs_flat)[I_intr]
end
@printf("Mean absolute error: %1.1f m\n", mean(abs.(dif)))
err_mean = mean(abs.(dif[I_no_ocean])*100 ./ obs_flat)
@printf("Mean abs error relative to true elevation: %1.3f %%\n", err_mean)

# retrieve matrix of reconstructed DEM
dem_rec             = zeros(nx*ny)
dem_rec[I_no_ocean] = x_rec
dem_rec_mat         = reshape(dem_rec,nx,ny)

# set values at < 400 m elevation (ixx) equal to aerodem
dem_rec_mat[ixx]   .= obs_orig[ixx]
# set very small values to zero as they are most likely ocean
dem_rec_mat[dem_rec_mat .< 5.] .= 0.

# smooth out (lots of jumps at data gaps where aerodem is enforced in adjacent pixels)
dem_smooth = mapwindow(median, dem_rec_mat, (5,5))

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
