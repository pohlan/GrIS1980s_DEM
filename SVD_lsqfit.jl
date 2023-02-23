using svd_IceSheetDEM, NetCDF, Statistics, LinearAlgebra, Glob, PyPlot, Printf, JLD2, TSVD
close("all")

# parameters and paths to adjust
r           = 40            # truncation of the SVD
λ           = 10^7          # regularization
model_files = glob("data/usurf_ex_gris_g1800m_v5_RAGIS_id_*.nc") # 60 time steps in total, 27 files

# load masks
I_no_ocean, I_data, I_intr = load("data/aerodem_data_indices.jld2","I_not_ocean","I_marg","I_intr")

# load model data
Data_all, nx, ny = read_model_data(which_files=2:27,tsteps=2:60;model_files)
# centering
Data       = Data_all[I_no_ocean, :]   # remove cells where there is ocean, saves half of the space
Data_mean  = mean(Data, dims=2)
Data_centr = Data .- Data_mean

# load observations
obs        = ncread(model_files[1],"usurf")[:,:,1]            # reconstruct a model geometry
# obs       = ncread("data/aerodem_gris_g1800.nc", "Band1")   # actual aerodem data
obs_flat   = reshape(obs, nx * ny, 1)[I_no_ocean]
x_data     = obs_flat .- Data_mean

# compute SVD
println("Computing the SVD..")
U, Σ, V = tsvd(Data_centr, r)

# solve the lsqfit problem
println("Solving the least squares problem..")
UΣ            = U*diagm(Σ)
U_A, Σ_A, V_A = svd(UΣ[I_data,:])
D             = transpose(diagm(Σ_A))*diagm(Σ_A) + λ*I
v_rec         = V_A * D^(-1) * transpose(diagm(Σ_A)) * transpose(U_A) * x_data[I_data]
x_rec         = U*diagm(Σ)*v_rec .+ Data_mean

dif              = zeros(nx*ny)
dif[I_no_ocean] .= x_rec .- obs_flat
# dif[I_intr] .= 0.
@printf("L2 error: %1.1f\n", norm(dif))

figure()
p = pcolormesh(reshape(dif,nx,ny)',cmap="bwr"); colorbar(label="[m]"); p.set_clim(-10,10)
title("reconstructed - true")
