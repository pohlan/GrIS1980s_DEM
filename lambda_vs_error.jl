include("SVD_lsqfit.jl")

errs = []
λs   = 10. .^(3:0.5:6.5)
for λ in λs
    D             = transpose(diagm(Σ_A))*diagm(Σ_A) + λ*I
    v_rec         = V_A * D^(-1) * transpose(diagm(Σ_A)) * transpose(U_A) * x_data[I_data]
    x_rec         = U*diagm(Σ)*v_rec .+ Data_mean
    # calculate error
    dif                      = zeros(nx*ny)
    dif[I_no_ocean[I_data]] .= (x_rec .- obs_flat)[I_data]
    err_mean = mean(abs.(dif))
    push!(errs, err_mean)
end

figure()
plot(λs, errs)
xscale("log")
xlabel("λ")
ylabel("mean absolute error")
title("resolution = $res m")
savefig("output/lambda_vs_error.jpg")
