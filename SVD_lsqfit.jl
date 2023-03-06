using svd_IceSheetDEM, NetCDF, Statistics, LinearAlgebra, Glob, PyPlot, Printf, JLD2, TSVD

# retrieve command line arguments
parsed_args = parse_commandline(ARGS)
r        = parsed_args["r"]        # truncation of the SVD
λ        = parsed_args["λ"]        # regularization
res      = parsed_args["res"]      # resolution
filepath = parsed_args["filepath"]
  
F        = Float32 # Julia default is Float64 but that kills the process for the full training data set if r is too large

# load masks
I_no_ocean, I_data, I_intr = load("data/indices_aerodem_g" * res * "m.jld2","I_not_ocean","I_marg","I_intr")

# load model data
model_files = glob(joinpath(filepath,"usurf_ex_gris_g" * res * "*_id_*YM.nc"))
Data_all, nx, ny = read_model_data(;F,model_files)
# centering
Data       = Data_all[I_no_ocean, :]   # remove cells where there is ocean, saves half of the space
Data_mean  = mean(Data, dims=2)
Data_centr = Data .- Data_mean

# compute SVD
println("Computing the SVD..")
if r === nothing
    U, Σ, V = svd(Data_centr)      # full SVD
else
    U, Σ, V = tsvd(Data_centr, r)
end
# load observations
# obs        = ncread(model_files[1],"usurf")[:,:,1]            # reconstruct a model geometry
obs_file = filepath * "aerodem_g" * res * "m_geoid_corrected_1978_1987_mean.nc"; obs = ncread(obs_file, "surface_altitude")
# obs_file = filepath * "bedmachine_" * res *  ".nc"; obs = ncread(obs_file, "surface_altitude")

obs_flat   = F.(reshape(obs, nx * ny, 1)[I_no_ocean])
x_data     = obs_flat .- Data_mean

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
    filename = "output/rec_" * obs_file[6:13] * "_r_$r"*"_lambda_$λ"*"_g$res.nc"
    varname  = "usurf"
    attribs  = Dict("units"   => "m",
                    "data_min" => 0.0)
    x    = ncread(model_files[1], "x")
    y    = ncread(model_files[1], "y")
    time = ncread(model_files[1], "time")
    data_rec = zeros(nx*ny)
    data_rec[I_no_ocean] = x_rec
    if isfile(filename)
        run(`rm $filename`)
    end
    nccreate(filename, varname, "x", x, "y", y, "time", 1, atts=attribs)
    ncwrite(reshape(data_rec,nx,ny,1), filename, varname)

    figure()
    p = pcolormesh(reshape(dif,nx,ny)',cmap="bwr"); colorbar(label="[m]"); p.set_clim(-50,50)
    title("reconstructed - true")
    savefig(filename[1:end-3]*".jpg")
end
