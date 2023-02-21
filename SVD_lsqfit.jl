using NetCDF, Statistics, LinearAlgebra, Glob, PyPlot, Printf, JLD2, TSVD
close("all")
# load model data
include("read_in.jl")
model_files = glob("data/usurf_ex_gris_g1800m_v5_RAGIS_id_*.nc") # 60 time steps in total, 27 files
Data, x, y, t, nx, ny = read_model_data(files=1:10,tsteps=1:60;model_files)
Data = permutedims(Data)

# center the data
Data_mean = mean(Data, dims=2)
Data_centr = Data .- Data_mean

println("Computing the SVD..")
U, S, V = tsvd(Data_centr, 60)

# load observations
obs       = ncread(model_files[27],"usurf")[:,:,1]
# obs       = ncread("data/aerodem_gris_g1800.nc", "Band1")
obs_flat  = reshape(obs, nx * ny, 1)
b         = obs_flat .- Data_mean
Id, I_intr = load("data/data1980s_indices.jld2","I_marg","I_intr")

# solve the lsqfit problem
λ     = 10^7
denom = S.^2 .+ λ
US    = U*diagm(S)
w     = transpose(US[Id,:]) * b[Id] ./ denom
# U_A, Σ_A, V_A = svd(US[Id,:])
# w     = V_A * diagm(Σ_A) * transpose(U_A) * b[Id] ./ denom  # this gives the exact same result as the first option to calculate w
x̂     = U*diagm(S)*w .+ Data_mean

dif   = x̂ .- obs_flat
# dif[I_intr] .= 0.
@printf("L2 error: %1.1f\n", norm(dif))

figure()
p = pcolormesh(reshape(dif,nx,ny)',cmap="bwr"); colorbar(label="[m]"); p.set_clim(-10,10)
title("reconstructed - true")
