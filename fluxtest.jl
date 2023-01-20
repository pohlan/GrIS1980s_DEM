## SIMPLE LINE FIT ##
# (example from https://fluxml.ai/Flux.jl/stable/models/overview/ but modified slightly since train! seems to want a params object as a second argument not the whole model)

# Dense(6 => 6) more efficient than Dense(1 => 1) with transposed x vector, both in terms of iterations and time
# For larger n_obs, LayerNorm seems to help! otherwise the loss becomes a NaN (also NaN with BatchNorm)
# not sure yet whether the additional hidden layers change much

using Flux, Statistics

actual(x) = 4x^2 + 2

n_obs = 200
n_train = Int(n_obs / 2)
x_train, x_test = hcat(1:n_train), hcat(n_train+1:n_obs)
y_train, y_test = actual.(x_train), actual.(x_test)

n_mid   = 50
m1 = Dense(n_mid=>n_mid)
predict = Chain(Dense(n_train=>n_mid),  LayerNorm(n_mid  , relu),
                Dense(n_mid=>n_mid),    LayerNorm(n_mid  ), Dropout(0.5),      # x -> m1(x) .+ Dense(n_mid => n_mid)(x),
                Dense(n_mid=>n_mid),    LayerNorm(n_mid  ), Dropout(0.1),
                # Dense(n_mid=>n_mid),    LayerNorm(n_mid  ),
                # Dense(n_mid=>n_mid),    LayerNorm(n_mid  ),
                # Dense(n_mid=>n_mid),    LayerNorm(n_mid  ),
               Dense(n_mid => n_train), LayerNorm(n_train),
                )

loss(x, y) = mean(abs2.(predict(x) .- y))
opt = Descent(0.1)
data = [(x_train, y_train)]
ps   = Flux.params(predict)

loss_0 = loss(x_train, y_train)
println("initial loss: $loss_0")
n_iters = 100000
tic = Base.time()
for k = 1:n_iters
    Flux.train!(loss, ps, data, opt)
end
toc = Base.time() - tic
loss_k = loss(x_train, y_train)
loss_test = loss(x_test, y_test)
println("Loss at iter=$n_iters: $loss_k after $toc seconds.")
println("Loss of test set: $loss_test.")
