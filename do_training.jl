function flux_training(;x_train,y_train,SV,λ,n_epochs=100)
    q = size(SV,1)
    model = Dense(q, 1)
    ps = Flux.params(model)
    penalty() = λ * mean(abs2, model.weight)
    loss(x, y) = Flux.Losses.mse(model(x), y) + penalty()
    opt = Adam(0.1)

    loss_0 = loss(x_train, y_train)
    println("Initial loss: $loss_0")

    data = [(x_train, y_train)]

    train_loader = Flux.DataLoader((x_train, y_train), batchsize = 128, shuffle = true)
    Flux.train!(loss, ps, ncycle(train_loader, n_epochs), opt)

    loss_e = loss(x_train, y_train)
    println("Final loss: $loss_e")
    M = model.weight * SV


    # Ms = copy(reshape(M, nx, ny, 1))
    # filename = "dem_ragis_q_$(q)_K_$(K).nc"
    # varname  = "usurf"
    # attribs  = Dict("units"   => "m",
    #                 "data_min" => 0.0)
    # nccreate(filename, varname, "x", x, "y", y, "time", 1, atts=attribs)
    # ncwrite(Ms, filename, varname)
    return M,loss_e
end


function optim_training(;x_train,y_train,US,λ)
    q = size(US,2)
    loss(V) = norm(x_train*V .- y_train) + λ * norm(V)
    V0 = zeros(q,1)
    res = optimize(loss, V0, BFGS())
    Vhat = Optim.minimizer(res)
    total_loss   = loss(Vhat)
    L2_loss = norm(x_train*Vhat .- y_train)
    M = US*Vhat
    return M, total_loss, L2_loss
end
