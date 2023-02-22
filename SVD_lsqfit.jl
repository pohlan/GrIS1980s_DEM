using NetCDF, Statistics, LinearAlgebra, Glob, PyPlot, Printf, JLD2, TSVD
close("all")

# parameters
r = 80            # truncation of the SVD
λ = 10^7          # regularization

# load model data
include("read_in.jl")
model_files = glob("data/usurf_ex_gris_g1800m_v5_RAGIS_id_*.nc") # 60 time steps in total, 27 files
Data, x, y, t, nx, ny = read_model_data(files=2:27,tsteps=2:60;model_files)
Data = permutedims(Data)

# center the data
Data_mean = mean(Data, dims=2)
Data_centr = Data .- Data_mean

println("Computing the SVD..")
U, Σ, V = tsvd(Data_centr, r)

# load observations
obs        = ncread(model_files[1],"usurf")[:,:,1]            # reconstruct a model geometry
# obs       = ncread("data/aerodem_gris_g1800.nc", "Band1")   # actual aerodem data
obs_flat   = reshape(obs, nx * ny, 1)
x_data     = obs_flat .- Data_mean
if !isfile("data/aerodem_data_indices.jld2")
    include("save_indices.jl")
end
Id, I_intr = load("data/data1980s_indices.jld2","I_marg","I_intr")

# solve the lsqfit problem
UΣ            = U*diagm(Σ)
U_A, Σ_A, V_A = svd(UΣ[Id,:])
D             = transpose(diagm(Σ_A))*diagm(Σ_A) + λ*I
v_rec         = V_A * D^(-1) * transpose(diagm(Σ_A)) * transpose(U_A) * x_data[Id] # this gives the exact same result as the first option to calculate w
x_rec         = U*diagm(Σ)*v_rec .+ Data_mean

dif   = x_rec .- obs_flat
# dif[I_intr] .= 0.
@printf("L2 error: %1.1f\n", norm(dif))

figure()
p = pcolormesh(reshape(dif,nx,ny)',cmap="bwr"); colorbar(label="[m]"); p.set_clim(-10,10)
title("reconstructed - true")
