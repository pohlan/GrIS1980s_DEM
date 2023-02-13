using Flux, TSVD, IterTools, Statistics, LinearAlgebra,
      Glob, NetCDF, PyPlot, Printf

# reading in model data
println("Reading in files...")
model_files = glob("data/usurf_ex_gris_g1800m_v5_RAGIS_id_*.nc") # 60 time steps in total, 27 files
t1, tend = 1, 30
nf, nt   = 18, 2  # number of files and time steps
d1   = ncread(model_files[1],"usurf"); nx, ny = size(d1[:,:,1])
Data = zeros(2*nf, nx*ny)
for (k, file) in enumerate(model_files[1:nf])
    d     = ncread(file, "usurf")[:,:,[t1, tend]]
    dflat = permutedims(reshape(d, nx*ny, nt))
    Data[(k-1)*2+1:2k,:] = dflat
end

# compute SVD
println("Computing SVD...")
Data_mean  = mean(Data,dims=1)
Data_centr = Data .- Data_mean
q = 2
U, S, V = tsvd(Data_centr, q)

# training and validation data
nf_train = 16   # number of files used for training
x_train, x_val  = permutedims(U[1:nt:nt*nf_train, :]), permutedims(U[nf_train*nt+1:nt:nf*nt,  :])
y_train, y_val  = permutedims(U[nt:nt:nt*nf_train,:]), permutedims(U[(nf_train+1)*nt:nt:nf*nt,:])

# tuning parameters
λ         = 7
lr        = 1e-2
n_epochs  = 100

# prepare for training
n_mid     = 20
model     = Chain(Dense(q => n_mid, leakyrelu),
                  Dense(n_mid => n_mid, leakyrelu),
                  Dense(n_mid => q))
ps        = Flux.params(model)
loss(x,y) = norm(model(x)-y,1) + λ * norm(ps)
opt       = Adam(lr)
data      = [(x_train, y_train)]

loss_0 = loss(x_train, y_train)
println("Initial loss: $loss_0.")

# do the training
train_loader = Flux.DataLoader((x_train, y_train), batchsize = 10, shuffle = true)
Flux.train!(loss, ps, ncycle(train_loader, n_epochs), opt)

loss_end = loss(x_train, y_train)
loss_val = loss(x_val,   y_val  )
println("Final loss: $loss_end.")
println("Loss val set: $loss_val.")
