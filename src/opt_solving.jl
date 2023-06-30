using Printf, Statistics, LinearAlgebra, TSVD

function solve_lsqfit(Data_all::Matrix{T}, obs::Matrix{T}, λ::T, r::Int, I_no_ocean::Vector{Int64}, I_obs::Vector{Int64}) where T<:Real
    nx, ny = size(obs)
    # centering model data
    Data       = Data_all[I_no_ocean, :]  # remove cells where there is ocean, saves half of the space
    Data_mean  = mean(Data, dims=2)
    Data_centr = Data .- Data_mean
    # centering observations with model mean
    obs_flat   = obs[I_no_ocean]
    x_data     = obs_flat .- Data_mean

    # compute SVD
    println("Computing the SVD..")
    U, Σ, _ = svd(Data_centr)

    # solve the lsqfit problem
    println("Solving the least squares problem..")
    UΣ            = U*diagm(Σ)
    U_A, Σ_A, V_A = svd(UΣ[I_obs,:])
    D             = transpose(diagm(Σ_A))*diagm(Σ_A) + λ*I
    v_rec         = V_A * D^(-1) * transpose(diagm(Σ_A)) * transpose(U_A) * x_data[I_obs]
    x_rec         = U*diagm(Σ)*v_rec .+ Data_mean

    # calculate error and print
    dif                     = zeros(T, nx,ny)
    dif[I_no_ocean[I_obs]] .= x_rec[I_obs] .- obs_flat[I_obs]
    err_mean         = mean(abs.(dif[I_no_ocean[I_obs]]))
    @printf("Mean absolute error: %1.1f m\n", err_mean)

    # retrieve matrix of reconstructed DEM
    dem_rec             = zeros(T, nx,ny)
    dem_rec[I_no_ocean] = x_rec

    return dem_rec, dif, err_mean
end
