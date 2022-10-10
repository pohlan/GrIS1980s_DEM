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

Idx = vec(u1 .> 0)

t = 150:13:292
x = Float32[]
x = vcat(x, vec(usurf[:,:,t[1]]))
for i in 2:length(t)
    x = hcat(x, vec(usurf[:,:,t[i]]))
end

x_mean = vec(mean(usurf,dims=3))
x_std  = vec(std(usurf,dims=3))
Xbar   = (x .- x_mean) ./ (x_std .+ eps())
R      = (vec(u1) .- x_mean) ./ (x_std .+ eps())
U, S, V = svd(Xbar);

# optimization problem
λ_u = 10.0    # σ^2 / σ_u^2
λ_v =  1.0    # σ^2 / σ_v^2

    # descent algorithm ...
    function f(Vi_)
        # Predicted matrix
        G = U * Diagonal(S) * Vi_ # Vi: column vector of VT    #
        # Data misfit
        E_misfit = sum((R[Idx] .- G[Idx]).^2)           # mean?
        # Regularization
        # E_reg = λ_u*sum(U.^2) + λ_v*sum(V_.^2)   # divide by len(Rbar)?
        # Sum to form total cost
        E = E_misfit #+ E_reg #  + E_space + E_time
        return E
    end
    return f, Xbar, u1, R, x_mean, x_std, testgap
end

f, Xbar, u1, R, x_mean, x_std, testgap = make_objective(usurf)

D = 6   # plot(S,yaxis=:log)
# U_ = 1e-2 .* randn(size(R,1),D)Infil
Vi_ = 1e-2 .* randn(size(Xbar,2),1)
f(Vi_)
out = optimize(f, Vi_)

v_est = out.minimizer
U, S, V = svd(Xbar)
xbar_est = U * Diagonal(S) * v_est
x_est = reshape(xbar_est .* x_std .+ x_mean,size(u1))
heatmap(((x_est.-usurf[:,:,1]).*isnan.(testgap))')
