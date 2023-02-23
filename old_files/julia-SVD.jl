using Flux, Optim, TSVD, IterTools, Statistics, LinearAlgebra,
      Glob, NetCDF, PyPlot, Printf, JLD2
close("all")

include("read_in.jl")
include("do_training.jl")

# load indices of data availability in 1980's dataset
Id, I_intr = load("data/data1980s_indices.jld2","I_marg","I_intr")

# load model data
model_files = glob("data/usurf_ex_gris_g1800m_v5_RAGIS_id_*.nc") # 60 time steps in total, 27 files
Data, x, y, t, nx, ny = read_model_data(files=1:15,tsteps=13:60;model_files)
# DEM to be reconstructed
t_y = 1
obs_artif = ncread(model_files[27],"usurf")[:,:,t_y]
R = reshape(obs_artif, 1, nx * ny)
# which_obs = 2
# obs_file = ["data/aerodem_1978_1987_wgs84_g1800m.nc",
#             "data/bamber_gris_g1800.nc"][which_obs]              # decide which one to use
# obs = ncread(obs_file, ["surface_altitude","Band1"][which_obs])
# R   = reshape(obs, 1, length(obs))

# center the data
Data_mean = mean(Data, dims=1)
Data_centr = Data .- Data_mean

q = 100
println("Computing svd...")
R_centr = R .- Data_mean         # center R as well
U, S, V = tsvd(Data_centr,q)
SV      = permutedims(V*diagm(S))
x_train = SV[:,Id]
_,n_obs = size(x_train)
y_train = reshape(R_centr[Id], 1, n_obs)

λs = [12]  
total_losses = zeros(length(λs))
errs_all     = zeros(length(λs))
errs_inn     = zeros(length(λs))
Ms       = [] 
Uhats    = []
###########################################
for (j,λ) in enumerate(λs)
    @printf("Training for λ = %1.1e, q = %d ...\n",λ,q)

    # do training
    tic = Base.time()
    # M, total_loss, Uhat = optim_training(;x_train,y_train,SV,λ)   # far too slow for bigger problems! (even with LBFGS)
    M, total_loss, Uhat = flux_training(;x_train,y_train,SV,λ,n_epochs=100)  
    toc = Base.time() - tic; println("  Elapsed time: $toc s.")
    M = M .+ Data_mean
    M[findall(x->x<0, M)] .= 0
    push!(Ms, M)
    push!(Uhats, Uhat)

    # save error
    total_losses[j] = total_loss
    errs_all[j] = norm(M.-R)
    errs_inn[j] = norm(M[I_intr] .- R[I_intr])
    @printf("  L2 norm error: %1.1f\n", errs_all[j])
end
###########################################

# plot difference between R (data) and Ms (fit)
dif = Ms[1] .- R
dif_rh = reshape(dif, nx,ny)
# dif_rh[obs.==-9999] .= NaN
dif_inner = zeros(size(R)); dif_inner[I_intr] .= dif[I_intr]
dif_outer = zeros(size(R)); dif_outer[Id] .= dif[Id]
figure()
p = pcolormesh(dif_rh',cmap="bwr"); colorbar(label="[m]"); p.set_clim(-10,10)
title("reconstructed - true")

# plot U's vs reconstructed Uhat 
figure()
for i in 1:27
    plot(U[i,:])
end
# plot(U[1,:])
plot(Uhats[end][:],lw=2, "black")
