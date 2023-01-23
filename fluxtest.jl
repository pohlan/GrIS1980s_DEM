using Flux, Statistics, IterTools, PyPlot

actual(x) = sqrt(x) # 4x^2 + 2

n_obs = 100
n_train = Int(n_obs / 2)
x_train, x_test = hcat(1:n_train...), hcat(n_train+1:n_obs...)
y_train, y_test = actual.(x_train), actual.(x_test)

n_mid   = 80
predict = Chain(Dense(1 => n_mid,     leakyrelu),
                Dense(n_mid => n_mid, leakyrelu),
                Dense(n_mid => n_mid, leakyrelu),
                Dense(n_mid=>1))
loss(x, y) = mean(abs2.(predict(x) .- y))
opt = Adam(1e-5)   # Descent(1e-5) not efficient, Adam much better!
data = [(x_train, y_train)]
train_loader = Flux.DataLoader((x_train, y_train), batchsize=10)
ps   = Flux.params(predict)

n_iters = 10000
loss_0 = loss(x_train, y_train)
println("initial loss: $loss_0")
tic = Base.time()
Flux.train!(loss,ps,ncycle(train_loader, n_iters), opt)
toc = Base.time() - tic
loss_end = loss(x_train, y_train)
loss_test = loss(x_test, y_test)
println("Loss at iter=$n_iters: $loss_end after $toc seconds.")
println("Loss of test set: $loss_test.")

# plot
close("all")
figure()
plot(x_train[:], predict(x_train)[:], label="training")
plot(x_test[:],  predict(x_test)[:],  label="validation")
plot([x_train[:];x_test[:]], [y_train[:]; y_test[:]], label="true")
legend()

# n_mid=20, batchsize=20, 2 hidden layers
# initial loss: 159.6740676990695
# Loss at iter=10000: 0.3666153805314159 after 8.119597911834717 seconds.
# Loss of test set: 14.285877127979756.

# n_mid=50, batchsize=20, 2 hidden layers
# initial loss: 354.0473716294248
# Loss at iter=10000: 0.03171823730434913 after 20.44066286087036 seconds.
# Loss of test set: 5.167972807210313

# n_mid=20, batchsize=20, 3 hidden layers
# initial loss: 87.56930224385452
# Loss at iter=10000: 0.06153866724325246 after 11.594810009002686 seconds.
# Loss of test set: 8.634141272600884

# n_mid=20, batchsize=5, 2 hidden layers
# initial loss: 192.97638417914575
# Loss at iter=10000: 0.001521786341840349 after 38.93712115287781 seconds.
# Loss of test set: 1.1089670876213196

# n_mid=100, batchsize=20, 2 hidden layers
# initial loss: 15.606674356185728
# Loss at iter=10000: 0.00013825186262780425 after 36.58544301986694 seconds.
# Loss of test set: 3.0275624997882447.

# n_mid=20, batchsize=20, 2 hidden layers
# initial loss: 12.047779718908217
# Loss at iter=50000: 0.0002722780206006659 after 51.98075008392334 seconds.
# Loss of test set: 0.555740149167875
