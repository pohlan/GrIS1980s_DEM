using NetCDF, Statistics, LinearAlgebra, Glob, PyPlot, Printf

# load model data
include("read_in.jl")
model_files = glob("data/usurf_ex_gris_g1800m_v5_RAGIS_id_*.nc") # 60 time steps in total, 27 files
Data, x, y, t, nx, ny = read_model_data(files=1:10,tsteps=1:60;model_files)
Data = permutedims(Data)

# center the data
Data_mean = mean(Data, dims=2)
Data_centr = Data .- Data_mean

println("Computing the SVD..")
U, S, V = svd(Data_centr)

# load observations
obs       = ncread(model_files[27],"usurf")[:,:,1]
obs_flat  = reshape(obs_artif, nx * ny, 1)
b         = obs_flat .- Data_mean

# solve the lsqfit problem
λ     = 0
denom = S.^2 .+ λ
w     = S .* transpose(U)*b ./ denom
x̂     = U*diagm(S)*w .+ Data_mean

dif   = x̂ .- obs_flat
@printf("L2 error: %1.1f\n", norm(dif))

figure()
p = pcolormesh(reshape(dif,nx,ny)',cmap="bwr"); colorbar(label="[m]"); p.set_clim(-10,10)
title("reconstructed - true")
