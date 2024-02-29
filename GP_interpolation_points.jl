using svd_IceSheetDEM
using NCDatasets, Interpolations, DataFrames, CSV, ProgressMeter, GeoStats, LsqFit, Distributed
addprocs(3)
@everywhere using StatsBase
import Plots

gr = 600

# masks
ds_mask = NCDataset("data/gris-imbie-1980/imbie_mask_g$(gr).nc")["Band1"]
bedm_mask = NCDataset("data/bedmachine/bedmachine_g$(gr).nc")["mask"]

# GrIMP
grimp = NCDataset("data/bedmachine/bedmachine_g$(gr).nc")["surface"]

# aero
ds_aero = NCDataset("data/grimp_minus_aero_g$(gr).nc")
dh_aero_all = ds_aero["dh"][:,:,1]
dh_aero_all = dh_aero_all[:,end:-1:1]
x  = sort(ds_aero["x"])
y  = sort(ds_aero["y"])

x_i    = 1:length(x)
y_i    = 1:length(y)
# x_i    = 750:1000
# y_i    = 150:1300
nx, ny = length(x_i), length(y_i)
idx_aero = findall(.!ismissing.(dh_aero_all[x_i,y_i]) .&& .!ismissing.(ds_mask[x_i,y_i]) .&& (bedm_mask[x_i,y_i] .!= 1))

df_aero_box = DataFrame(:x       => Float64.(repeat(x[x_i], 1, ny)[idx_aero]),
                        :y       => Float64.(repeat(y[y_i]', nx, 1)[idx_aero]),
                        :h_grimp => Float64.(grimp[x_i,y_i][idx_aero]),
                        :dh      => Float64.(dh_aero_all[x_i,y_i][idx_aero]),
                        :idx     => idx_aero)

# atm
df_atm = CSV.read("atm_grimp_interp_pts.csv", DataFrame)

# take averages
# df_atm_avg = DataFrame(names(df_atm) .=> [Float64[] for t in eachcol(df_atm)])
# nd = 600
# iblob = 1:nd:size(df_atm,1)
# for i in iblob
#     iend = min(i+500,size(df_atm,1))
#     push!(df_atm_avg,mean.(eachcol(df_atm[i:iend,:])))
# end


x_min, x_max = extrema(x[x_i])
y_min, y_max = extrema(y[y_i])
mindist = 5e4 # minimum distance between aerodem and atm points
function is_in_box(x, y)::Bool
    dist_to_aero = minimum(sqrt.((x.-df_aero_box.x).^2 .+ (y.-df_aero_box.y).^2))
    dist_to_aero = minimum(pairwise(Distances.euclidean, [x y], [df_aero_box.x df_aero_box.y], dims=1)[:])
    # is_in_x = x_min < x < x_max
    # is_in_y = y_min < y < y_max
    is_above_mindist = dist_to_aero > mindist
    return is_above_mindist
end
id = unique(sort(rand(1:size(df_atm,1), 20000)))
println("selecting flightline values...")
df_atm_box = filter([:x,:y] => is_in_box, df_atm[id,:])

# remove outliers
atm_to_delete = findall(abs.(df_atm_box.dh) .> 3 .* mad(df_atm_box.dh))
deleteat!(df_atm_box, atm_to_delete)

aero_to_delete = findall(abs.(df_aero_box.dh) .> 3 .* mad(df_aero_box.dh) .|| df_aero_box.h_grimp .< 50)
deleteat!(df_aero_box, aero_to_delete)

# plot again without outliers
dh_plot = Array{Union{Float32,Missing}}(missing, (nx,ny))
dh_plot[df_aero_box.idx] = df_aero_box.dh
Plots.heatmap(x[x_i],y[y_i], dh_plot', cmap=:bwr, clims=(-10,10))
Plots.scatter!(df_atm_box.x, df_atm_box.y, marker_z=df_atm_box.dh, color=:bwr, markersize=2, markerstrokewidth=0, legend=false)
Plots.savefig("dplot.png")

df_all = vcat(df_aero_box, df_atm_box, cols=:intersect)

# remove mean trend as a fct of elevation
dh_binned, ELV_bin_centers = svd_IceSheetDEM.bin_equal_sample_size(df_all.h_grimp, df_all.dh, 16000)  # 16000
itp_bias = linear_interpolation(ELV_bin_centers, mean.(dh_binned),  extrapolation_bc = Interpolations.Flat())
Plots.scatter(ELV_bin_centers, mean.(dh_binned))
ELV_plot = minimum(df_all.h_grimp):10:maximum(df_all.h_grimp)
Plots.plot!(ELV_plot, itp_bias(ELV_plot), label="linear interpolation", legend=:bottomright)
Plots.savefig("c1plot.png")

# standardize
itp_std = linear_interpolation(ELV_bin_centers, std.(dh_binned),  extrapolation_bc = Interpolations.Flat())
Plots.scatter(ELV_bin_centers, std.(dh_binned))
ELV_plot = minimum(df_all.h_grimp):10:maximum(df_all.h_grimp)
Plots.plot!(ELV_plot, itp_std(ELV_plot), label="linear interpolation", legend=:bottomright)
Plots.savefig("c2plot.png")

df_all.dh_detrend = (df_all.dh .- itp_bias(df_all.h_grimp)) ./ itp_std(df_all.h_grimp)

# table = (; Z=df_all.dh_detrend)
table = (; Z=df_all.dh_detrend)
# coords = [(df_all.x[i], df_all.y[i]) for i in 1:length(df_all.x)]
coords = [(df_all.x[i], df_all.y[i]) for i in 1:length(df_all.x)]
data = georef(table,coords)

println("Variogram...")
gamma = EmpiricalVariogram(data, :Z; nlags=400,  maxlag=6e5)

############################################
# function SEIso_kernel(x, params)
#     γ1 = params[1]^2 .* (1 .- exp.(-x.^2 ./ (2*params[2]^2)))
#     # γ2 = params[3]^2 .* (1 .- exp.(-x.^2 ./ (2*params[4]^2)))
#     # return γ1 .+γ2
# end
# function RQIso_kernel(x, params)
#     ℓ, σ, α = params
#     γ1 = σ^2 * (-x.^2 ./ (2*α*ℓ^2)).^(-α)
#     return γ1
# end
# function Mat12Iso_kernel(x, params)
#     ℓ, σ = params
#     cov = σ^2 .* exp.(-x./ℓ)
#     γ1 = σ^2 .- cov
#     return γ1
# end
# function Mat52Iso_kernel(x, params)
#     ℓ, σ = params
#     cov = σ^2 * (1 .+ √5 .*x./ℓ+5*x.^2/(3ℓ^2)) .* exp.(-√5 .*x./ℓ)
#     γ1 = σ^2 .- cov
#     return γ1
# end

# # ff = LsqFit.curve_fit(custom_kernel, gamma.abscissa, gamma.ordinate, [5, 2e4]) #, 5, 2e4]);
# Plots.scatter(gamma.abscissa, gamma.ordinate, xscale=:log10)
# Plots.plot!(gamma.abscissa, SEIso_kernel(gamma.abscissa, [1.0, 1e5]), label="SEIso")
# # Plots.plot!(gamma.abscissa, RQIso_kernel(gamma.abscissa, [1e4, 1.0, -0.5]))
# Plots.plot!(gamma.abscissa, Mat52Iso_kernel(gamma.abscissa, [4e4, 1.0]), label="Matern52")
# Plots.plot!(gamma.abscissa, Mat12Iso_kernel(gamma.abscissa, [3e4, 1.1]), label="Matern12")
##############################################




# GeoStats.interpolate

##################################


function custom_var(x, params)
    γ1 = ExponentialVariogram(range=params[1], sill=params[2], nugget=0.0)
    # γ2 = ExponentialVariogram(range=params[3], sill=params[4], nugget=0.0)
    f = γ1.(x) #.+ γ2.(x)
    return f
end
ff = LsqFit.curve_fit(custom_var, gamma.abscissa, gamma.ordinate, [5e4, 1.0]);

# julia> ff.param
# 2-element Vector{Float64}:
#  59344.64967915991
#      0.94160717609134
varg = ExponentialVariogram(range=ff.param[1], sill=ff.param[2]) # +
    #    ExponentialVariogram(range=ff.param[3], sill=ff.param[4])
Plots.scatter(gamma.abscissa, gamma.ordinate)
Plots.plot!(gamma.abscissa,varg.(gamma.abscissa), xscale=:log)


idx = findall(.!ismissing.(ds_mask[x_i,y_i]) .&& (bedm_mask[x_i,y_i] .!= 1))

df_idx = DataFrame(:x       => Float64.(repeat(x[x_i], 1, ny)[idx]),
                   :y       => Float64.(repeat(y[y_i]', nx, 1)[idx]),
                   :h_grimp => Float64.(grimp[x_i,y_i][idx]),
                   :idx     => idx)

# xi_predict = 700:1000; xip_min, xip_max = extrema(x[xi_predict])
# yi_predict = 250:600; yip_min, yip_max = extrema(y[yi_predict])
# ipred = findall(xip_min.<df_all.x .< xip_max .&& yip_min.<df_all.y .< yip_max)
# table_pred = (; Z=df_all.dh_detrend[ipred])

# # x0, xend = extrema(x[])
# # y0, yend = extrema(y[])
# grid = CartesianGrid((xip_min,yip_min), (xip_max,yip_max), dims=(length(xi_predict),length(yi_predict)))

table    = (; Z=Float32.(df_all.dh_detrend))
coords   = [(xi,yi) for (xi,yi) in zip(df_all.x, df_all.y)]
geotable = georef(table,coords)
grid     = CartesianGrid((x[1]-gr/2, y[1]-gr/2), (x[end]+gr/2, y[end]+gr/2), dims=(length(x), length(y)))

# solver = Kriging(:Z => (variogram=varg,))
# problem = EstimationProblem(data_pred, grid, :Z)
# solution = solve(problem, solver)
# mb = reshape(solution.Z, length(xi_predict),length(yi_predict))

model = Kriging(varg)
println("Kriging...")
tic = Base.time()
interp = geotable |> Interpolate(grid, model)
toc = Base.time() - tic
# tt = toc / 60
println(toc)
mb = reshape(interp.Z, length(xi_predict),length(yi_predict))

# Plots.heatmap(x[xi_predict], y[yi_predict], mb', cmap=:bwr, clims=(-10,10))

h_predict = vec(grimp[xi_predict, yi_predict]) .- (vec(mb).*itp_std(vec(grimp[xi_predict, yi_predict])) .+ itp_bias(vec(grimp[xi_predict, yi_predict])))
Plots.heatmap(x[xi_predict], y[yi_predict], reshape(h_predict, nx,ny)')
Plots.savefig("h_predict.png")

h_aero = grimp[xi_predict,yi_predict] .- dh_aero_all[xi_predict,yi_predict]
dd = reshape(h_predict,nx,ny) .- h_aero
d = Array{Union{Float32,Missing}}(missing, (nx,ny))
d[idx_aero] = dd[idx_aero]
Plots.heatmap(x[xi_predict], y[yi_predict], d', cmap=:bwr, clims=(-15,15))
Plots.savefig("diff_to_aero.png")


# GP

#########################

# m = MeanZero()
# kern  = Mat12Iso(log(3e4), log(1.1)) #+
#         # Mat12Iso(log(3e4), log(1.1))

# println("GP")
# gp = GP([df_all.x df_all.y]',df_all.dh_detrend,m,kern, log(1.0))
# println("optimize...")
# optimize!(gp)

# n_pts = 3000
# idx_u = rand(1:length(x_all), n_pts)
# gp = GaussianProcesses.SoR([x_all y_all]', [x_all[idx_u] y_all[idx_u]]', dh_detrend, m, kern, log(200));
# optimize!(gp_SOR)

# # predict
# println("predict...")
# idx_predict = findall(.!ismissing.(ds_mask[x_i,y_i]) .&& (bedm_mask[x_i,y_i] .!= 1) .&& .!ismissing.(grimp[x_i,y_i]))
# x_predict   = repeat(x[x_i],  1, ny)[idx_predict]
# y_predict   = repeat(y[y_i]', nx, 1)[idx_predict]
# elv_predict =         grimp[x_i,y_i][idx_predict]
# # s_predict   =   grimp_slope[x_i,y_i][idx_predict]
# mus, sig    = predict_y(gp, [vec(x_predict) vec(y_predict)]') # vec(elv_predict) vec(s_predict)]')
# dh_pred     =  Array{Union{eltype(dh_aero),Missing}}(missing, (nx,ny))
# dh_pred[idx_predict] .= mus .+ itp(elv_predict)
# Plots.heatmap(x[x_i], y[y_i], dh_pred', cmap=:bwr, clims=(-30,30)); Plots.savefig("bplot.png")
