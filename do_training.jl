function flux_training(;x_train,y_train,SV,λ,n_epochs=100)
    q = size(SV,1)
    model = Dense(q, 1)
    ps = Flux.params(model)
    penalty() = λ * mean(abs2, model.weight)
    loss(x, y) = Flux.Losses.mse(model(x), y) + penalty()
    opt = Adam(0.1)

    loss_0 = loss(x_train, y_train)
    println("  Initial loss: $loss_0")

    data = [(x_train, y_train)]

    train_loader = Flux.DataLoader((x_train, y_train), batchsize = 128, shuffle = true)
    Flux.train!(loss, ps, ncycle(train_loader, n_epochs), opt)

    loss_e = loss(x_train, y_train)
    println("  Final loss: $loss_e")
    M = model.weight * SV


    # Ms = copy(reshape(M, nx, ny, 1))
    # filename = "dem_ragis_q_$(q)_K_$(K).nc"
    # varname  = "usurf"
    # attribs  = Dict("units"   => "m",
    #                 "data_min" => 0.0)
    # nccreate(filename, varname, "x", x, "y", y, "time", 1, atts=attribs)
    # ncwrite(Ms, filename, varname)
    return M,loss_e,model.weight
end

function optim_training(;x_train,y_train,SV,λ)
    q       = size(SV,1)
    loss(U) = norm(U*x_train .- y_train) + λ * norm(U)
    U0      = zeros(1,q)
    loss_0  = loss(U0); println("  Initial loss: $loss_0.")
    res     = optimize(loss, U0, BFGS())
    Uhat    = Optim.minimizer(res)
    loss_e  = loss(Uhat); println("  Final loss: $loss_e.")
    M       = Uhat*SV
    return M, loss_e, Uhat
end
