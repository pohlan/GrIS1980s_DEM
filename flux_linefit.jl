## SIMPLE LINE FIT ##
# (example from https://fluxml.ai/Flux.jl/stable/models/overview/ but modified slightly since train! seems to want a params object as a second argument not the whole model)

using Flux, Statistics

actual(x) = 4x + 2

x_train, x_test = hcat(1:5...), hcat(6:10...)
y_train, y_test = actual.(x_train), actual.(x_test)

predict = Dense(1 => 1)
ps   = Flux.params(predict)

loss(x, y) = mean(abs2.(predict(x) .- y))
opt = Descent(0.08)   # for the default value η=0.1, the loss increases!
                      # Adam(0.3) is a bit better, but Adam(0.1) is potentially worse than Descent(0.08)
data = [(x_train, y_train)]

loss_0 = loss(x_train, y_train)
println("initial loss: $loss_0")

n_iters = 1000
for k = 1:n_iters
    Flux.train!(loss, ps, data, opt)
end

loss_k = loss(x_train, y_train)
loss_test = loss(x_test, y_test)
println("Loss at iter=$n_iters: $loss_k.")
println("Loss of test set: $loss_test.")
