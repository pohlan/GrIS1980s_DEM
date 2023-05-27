using svd_IceSheetDEM, NetCDF, Statistics, LinearAlgebra, Glob, PyPlot, Printf, TSVD, ImageFiltering
import ArchGDAL as AG

ARGS = ["--save", "--lambda", "1e6", "--obs", "data/aerodem_filtered_g1200m_geoid.nc", "--r", "80",
        "--obs_band_name", "surface"]

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

# save as nc file
if parsed_args["save"]
    mkpath("output/")
    println("Saving file..")
    logλ = Int(round(log(10, λ)))
    filename = "output/rec_lambda_1e$logλ"*"_g$res"*"_r$r.nc"
    varname  = "usurf"
    data_rec = zeros(nx*ny)
    data_rec[I_no_ocean] = x_rec
    data_matrix = reshape(data_rec,nx,ny)
    data_matrix[ixx] .= obs_orig[ixx]
    data_matrix[data_matrix .< 0.] .= 0.
    data_matrix = data_matrix[:,end:-1:1]  # a bit of a hack; turn Greenland 'upside down' so that it is correct in the final file
    model_dataset = AG.read(obs_file)
    AG.create(
        filename,
        driver = AG.getdriver(model_dataset),
        width  = AG.width(model_dataset),
        height = AG.height(model_dataset),
        nbands = 1,
        dtype  = Float32
    ) do raster
        AG.write!(raster, data_matrix, 1)
        AG.setgeotransform!(raster, AG.getgeotransform(model_dataset))
        AG.setproj!(raster, AG.getproj(model_dataset))
    end
    fn_compressed = filename[1:end-3] * "_compressed.nc"
    run(`gdal_translate -co "COMPRESS=DEFLATE" $filename $fn_compressed`)

    figure(figsize=(14,16))
    p = pcolormesh(reshape(dif,nx,ny)',cmap="bwr"); colorbar(label="[m]"); p.set_clim(-200,200)
    title("reconstructed - true")
    savefig(filename[1:end-3]*".jpg")
end
