using Flux, Plots, TSVD, IterTools, Statistics, LinearAlgebra, Compat, Glob, NetCDF

include("read_in.jl")
include("do_training.jl")

Data, I, I_inv, R = make_training_data(rtrain=1:8,tsteps=1:1,t_y=40)

qs = 2:15
losses_e = zeros(length(qs))
dhs      = zeros(length(qs))
for (i,λ) in enumerate([0 1e2 1e3 1e4 1e5 1e6 1e7])
    println("Training for λ = $λ ...")
###########################################
# center the data and compute the SVD
q = 5
Data_mean = mean(Data, dims=1)
Data_centr = Data .- Data_mean
U, S, V = tsvd(Data_centr, q)
SV = transpose(V * diagm(S))

# prepare training
x_train = transpose(SV[repeat(I, 1, q)])
_, n_obs = size(x_train)
R_centr = R .- Data_mean         # center R as well
y_train = reshape(R_centr[I], 1, n_obs)
###########################################


# do training
#λ = 1e4
M, loss_e = train_data(;x_train,y_train,SV,λ)
M = M .+ Data_mean
M[findall(x->x<0, M)] .= 0

losses_e[i] = loss_e
dhs[i] = Flux.Losses.mse(M[I_inv],R[I_inv])

end

# lloss
# qqs
# dhh

# filename = "dem_R.nc"
# jldsave(filename;M)



# "true" final field at the margin -> blend between data and minimization result ?
# smoothing where the interior touches I?


