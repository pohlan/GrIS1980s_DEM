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

I = vec(u1 .> 0)

t = 100:50:200
x = Float32[]
x = vcat(x, vec(usurf[:,:,t[1]]))
for i in 2:length(t)
    x = hcat(x, vec(usurf[:,:,t[i]]))
    I = hcat(I, trues(length(usurf[:,:,1]),1))
end

x_mean = mean(x,dims=2)
x_std  = std(x,dims=2)
R = (x .- x_mean) ./ (x_std .+ eps())

# optimization problem
λ_u = 10.0    # σ^2 / σ_u^2
λ_v =  1.0    # σ^2 / σ_v^2

    # descent algorithm ...
    function f(U_VT_)
        U_ = U_VT_[1:size(R,1),:]
        V_ = U_VT_[size(R,1)+1:end,:]'
        # Predicted matrix
        G = U_ * V_
        # Data misfit
        E_misfit = sum((R[I] .- G[I]).^2)           # mean?
        # Regularization
        E_reg = λ_u*sum(U_.^2) + λ_v*sum(V_.^2)   # divide by len(Rbar)?
        # Sum to form total cost
        E = E_misfit + E_reg #  + E_space + E_time
        return E
    end
    return f, R
end

f, R = make_objective(usurf)

U, S, V = svd(R);
D = 6   # plot(S,yaxis=:log)
U_ = 1e-2 .* randn(size(R,1),D)
V_ = 1e-2 .* randn(D,size(R,2))
x0 = [U_; V_']
# optimize(f, x0)    # does not work
