using Flux, Optim, TSVD, IterTools, Statistics, LinearAlgebra,
      Compat, Glob, NetCDF, PyPlot, Printf, DelimitedFiles, JLD2

include("read_in.jl")
include("do_training.jl")

I, I_inv = load("data/data1980s_indices.jld2","I","I_inv")

model_files = glob("data/usurf_gris_g1800m_v5_RAGIS_id_*_1980-1-1_2020-1-1_YM.nc") # 40 time steps in total, 9 files
# x-data
Data, nx, ny = xdata_training(files=1:8,tsteps=5:30;model_files)
# y-data:
t_y = 10
obs_artif = ncread(model_files[end],"usurf")[:,:,t_y]
R = reshape(obs_artif, nx * ny, 1);

# center the data
Data_mean = mean(Data, dims=2)
Data_centr = Data .- Data_mean

qs = [8] #4:10
λs = [10 50] # 80 100 300 500]
total_losses = zeros(length(qs), length(λs))
dhs      = zeros(length(qs), length(λs))
Ms       = zeros(length(qs), length(λs),nx*ny)
###########################################
for (i,q) in enumerate(qs)
    # compute the SVD
    U, S, V = tsvd(Data_centr, q)
    US = U * diagm(S)
    # prepare training
    x_train = US[repeat(I, 1, q)]
    n_obs, _ = size(x_train)
    R_centr = R .- Data_mean         # center R as well
    y_train = reshape(R_centr[I], n_obs, 1)

    for (j,λ) in enumerate(λs)
        @printf("Training for λ = %1.1e, q = %d ...\n",λ,q)

        # do training
        M, total_loss, L2_loss = optim_training(;x_train,y_train,US,λ)
        M = M .+ Data_mean
        M[findall(x->x<0, M)] .= 0

        # save error
        total_losses[i,j] = total_loss
        dhs[i,j]          = norm(M[I_inv].-R[I_inv],2)
        # dhs[i,j]        = norm(M.-R,2)
        Ms[i,j,:] = M
    end
end
###########################################

# close()
# l1_inner = [norm(Ms[1,i,I_inv] .- R[I_inv],1) for i = 1:length(dhs[:])]
# l2_inner = [norm(Ms[1,i,I_inv] .- R[I_inv],2) for i = 1:length(dhs[:])]
# l1_all   = [norm(Ms[1,i,:] .- R,1) for i = 1:length(dhs[:])]
# l2_all   = [norm(Ms[1,i,:] .- R,2) for i = 1:length(dhs[:])]
# subplot(2,2,1)
# plot(l1_inner,"-o"); title("l1_inner")
# subplot(2,2,2)
# plot(l2_inner,"-o"); title("l2_inner")
# subplot(2,2,3)
# plot(l1_all,"-o"); title("l1_all")
# subplot(2,2,4)
# plot(l2_all,"-o"); title("l2_all")

# plot result of optimization (fitted ice thickness)
pcolor(reshape(Ms[1,1,:],nx,ny))

# plot difference between R (data) and Ms (fit)
dif = R .- Ms[1,1,:]
dif[dif.==0] .= NaN
dif_inner = copy(dif); dif_inner[:] .= NaN; dif_inner[I_inv] .= dif[I_inv]
dif_outer = copy(dif); dif_outer[:] .= NaN; dif_outer[I] .= dif[I]
figure()
pcolor(reshape(dif_outer,nx,ny)',cmap="bwr",clim=(-150,150))
colorbar()
