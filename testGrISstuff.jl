using NetCDF, Plots, LinearAlgebra, Statistics, Optim
file = "../Downloads/ex_gris_g1800m_v3a_rcp_26_id_0_2008_2300.nc"
# ncinfo(file)
usurf = ncread(file, "usurf") # surface elevation


function make_objective(usurf)

testgap = ones(size(usurf[:,:,1]))
testgap[350:600,600:1150] .= NaN
testgap[290:520,800:1240] .= NaN
testgap[300:400,200:800]  .= NaN
testgap[400:500,430:600]  .= NaN
u1 = usurf[:,:,1] .* testgap
heatmap(u1')

Idx = vec(u1 .> 0)

t = 150:50:292
x = Float32[]
x = vcat(x, vec(usurf[:,:,t[1]]))
for i in 2:length(t)
    x = hcat(x, vec(usurf[:,:,t[i]]))
    Idx = hcat(Idx, trues(length(usurf[:,:,1]),1))
end

x_mean = mean(x,dims=2)
x_std  = std(x,dims=2)
R = (x .- x_mean) ./ (x_std .+ eps())

U, S, V = svd(R);

# optimization problem
λ_u = 10.0    # σ^2 / σ_u^2
λ_v =  1.0    # σ^2 / σ_v^2

    # descent algorithm ...
    function f(V_)
        # Predicted matrix
        G = U * Diagonal(S) * V_'
        # Data misfit
        E_misfit = sqrt(sum((R[Idx] .- G[Idx]).^2)) / sum(Idx)           # mean?
        # Regularization
        # E_reg = λ_u*sum(U.^2) + λ_v*sum(V_.^2)   # divide by len(Rbar)?
        # Sum to form total cost
        E = E_misfit #+ E_reg #  + E_space + E_time
        return E
    end
    return f, R
end

f, R = make_objective(usurf)

D = 6   # plot(S,yaxis=:log)
# U_ = 1e-2 .* randn(size(R,1),D)Infil
V_ = 1e-2 .* randn(size(R,2),size(R,2))
f(V_)
out = optimize(f, V_)
