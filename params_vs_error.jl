using svd_IceSheetDEM, NetCDF, PyPlot

# three parameters evaluated against mean absolute error:
# 1) lambda (regularization) vs error
# 2) number of ensembles vs error
# 3) r (SVD truncation) vs erro
# note that process is easily killed if all of the three above are done at once

const F = Float32

parsed_args = parse_commandline(ARGS)
training_data_files = parsed_args["training_data"]
imbie_path          = "data/gris-imbie-1980/"
aerodem_path        = "data/aerodem/"

x = ncread(training_data_files[1], "x")
const gr = Int(x[2] - x[1])   # assumes same grid size in both x and y direction

obs_file            = aerodem_path*"aerodem_rm-filtered_geoid-corr_g$(gr).nc"
imbie_mask          = imbie_path * "imbie_mask_g$(gr).nc"

Data_all, obs, I_no_ocean, I_obs, nx, ny, _ = svd_IceSheetDEM.prepare_problem(obs_file, imbie_mask, training_data_files, F)

###################
# lambda vs error #
###################

r = parsed_args["r"]
errs_λ = []
λs     = 10. .^(3:0.5:6.5)
for λ in λs
    _, err_mean, _ = svd_IceSheetDEM.solve_problem(Data_all, obs, I_no_ocean, I_obs, nx, ny, r, λ, F)
    push!(errs_λ, err_mean)
end

figure()
plot(λs, errs_λ)
xscale("log")
xlabel("λ")
ylabel("mean absolute error")
if r > 700
    title("resolution = $gr m, full SVD")
else
    title("resolution = $gr m, r = $r")
end
savefig("output/lambda_vs_error.jpg")


################################
# number of ensembles vs error #
################################

# nf     = length(training_data_files)
# nt     = Int(size(Data_all,2) / nf)
# λ      = parsed_args["λ"]
# r = parsed_args["r"]
# errs_m = []
# ms     = 5:2:nf
# for m in ms
#     Data = Data_all[:, 1:m*nt]
#     _, err_mean, _ = svd_IceSheetDEM.solve_problem(Data, obs, I_no_ocean, I_obs, nx, ny, r, λ, F)
#     push!(errs_m, err_mean)
# end

# figure()
# plot(ms, errs_m)
# # xscale("log")
# xlabel("# ensembles")
# ylabel("mean absolute error")
# logλ = Int(log10(λ))
# if r > 700
#     title("resolution = $gr m, λ = 1e$logλ, full SVD")
# else
#     title("resolution = $gr m, λ = 1e$logλ, r = $r")
# end
# savefig("output/n_ensembles_vs_error.jpg")


#########################
# truncation r vs error #
#########################

# λ      = parsed_args["λ"]
# errs_r = []
# rs     = [10, 50, 100, 200, 300, 400, 500, 800]
# for r in rs
#     _, err_mean, _ = svd_IceSheetDEM.solve_problem(Data_all, obs, I_no_ocean, I_obs, nx, ny, r, λ, F)
#     push!(errs_r, err_mean)
# end

# figure()
# plot(rs, errs_r)
# xscale("log")
# xlabel("SVD truncation r")
# ylabel("mean absolute error")
# logλ = Int(log10(λ))
# title("resolution = $gr m, λ = 1e$logλ")
# savefig("output/r_vs_error.jpg")
