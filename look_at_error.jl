using svd_IceSheetDEM, NetCDF, Interpolations, GeoStats, GaussianRandomFields, StatsBase
import Plots, StatsPlots
import LsqFit as LF



function bin_equal_sample_size(x, y, n_samples)
    nx = length(x)
    n_bins = ceil(Int, nx / n_samples)
    if nx % n_samples < 100
        n_bins -= 1
    end
    p = sortperm(x)

    out = Vector{Real}[]
    bin_centers = Real[]
    for b in 1:n_bins
        if b == n_bins
            i = (b-1)*n_samples+1:nx
        else
            i = (b-1)*n_samples+1:b*n_samples
        end
        push!(bin_centers, 0.5*(sort(x)[i[1]]+sort(x)[i[end]]))
        push!(out, y[p[i]])
    end
    return out, bin_centers
end

# load observations

obs_aero_file = "data/aerodem/aerodem_rm-filtered_geoid-corr_g600.nc"
rec_file      = "output/bedmachine1980_reconstructed_g600.nc"

rec = ncread(rec_file, "surface")
x   = ncread("output/bedmachine1980_reconstructed_g600.nc", "x")
y   = ncread("output/bedmachine1980_reconstructed_g600.nc", "y")
mask = ncread("output/bedmachine1980_reconstructed_g600.nc", "mask")
obs_aero = ncread(obs_aero_file, "Band1")
obs_ATM = ncread("data/ATM/ATM_geoid_corrected_g600.nc", "surface")


# calculate difference between reconstruction and observations (where available)
idx_0 = findall(vec(rec .!= 0) .&& vec(obs_aero .!= 0) .&& vec(mask .!= 1.0) .&& vec(rec .!= obs_aero))
idx_1 = findall(vec(rec .!= 0) .&& vec(obs_ATM .> 3000) .&& vec(mask .!= 1.0) .&& vec(obs_aero .== 0))

# obs = obs_aero

dh_0 = rec[idx_0] .- obs_aero[idx_0]
dh_1 = rec[idx_1] .- obs_ATM[idx_1]

obs = obs_aero
idx = idx_0
dh  = dh_0

# obs = vcat(vec(obs_aero), vec(obs_ATM))
# dh = vcat(dh_0, dh_1)
# idx = vcat(idx_0, idx_1 .+ length(obs_aero))

# std_dh = mad(dh_all[idx])
# filter!(x -> abs(dh_all[x]) < 10.0*std_dh, idx)
# dh = dh_all[idx]

# # calculate slope
# dfx = diff(diff(obs, dims=1), dims=1)
# dfy = diff(diff(obs, dims=2), dims=2)
# df = max.(abs.(dfx[:,2:end-1]), abs.(dfy[2:end-1,:]))

# seperate elevation differences into bins along slope of the corresponding observation
# dh_binned_df, df_bin_centers = get_bins(1000, dh, df[idx])
# nmad_df = std.(dh_binned_df)
# bias_df = median.(dh_binned_df)
# Plots.scatter(df_bin_centers, nmad_df, xlabel="slope", ylabel="nmad", legend=:false)
# Plots.scatter(df_bin_centers, bias_df, xlabel="slope", ylabel="bias [m]", legend=:false)
# StatsPlots.violin(dh_binned_df, linewidth=0.5, color=:darkblue, legend=:false, fillalpha=0.6, ylabel="elevation difference [m]")


##### ANALYSIS W.R.T. ELEVATION ######
### seperate elevation differences into bins along elevation of the corresponding observation
nsamples = 1000
dh_binned, bin_centers = bin_equal_sample_size(obs[idx], dh, nsamples)
idx_binned, _ = bin_equal_sample_size(obs[idx], idx, nsamples)

### filter per bin
for (i_bin, dh_bin) in zip(idx_binned, dh_binned)
    p = sortperm(abs.(dh_bin))
    n = length(dh_bin)
    nmad_ = StatsBase.mad(dh_bin)
    filter!(x -> abs(x) < 7.0*nmad_, dh_bin)
    deleteat!(i_bin, sort(p[length(dh_bin)+1:n]))   # get the indexes of the dh_bin in ascending order of their absolute value, and take the first n
end
# update dh and idx
dh = vcat(dh_binned...)
idx = vcat(idx_binned...)


### violin plot
StatsPlots.violin(dh_binned, linewidth=0.5, color=:darkblue, legend=:false, fillalpha=0.6, ylabel="elevation difference [m]")


### sigma
# calculate normalized mean absolute deviation (NMAD) of these bins
nmad = StatsBase.mad.(dh_binned)

## model the empirical dispersion as a function of elevation
# linear interpolation
lin_interp = linear_interpolation(bin_centers, nmad, extrapolation_bc=Interpolations.Line())
# curve fitting
@. linear_model(x,p) = p[1] + x*p[2]
fit_lin = LF.curve_fit(linear_model, bin_centers, nmad, [-4.5,0.015])
@. exp_model(x,p) = p[1]*exp(x*p[2])
fit = LF.curve_fit(exp_model, bin_centers, nmad, [15.,7e-4])
#plot
x_obs = minimum(obs[idx]):10:maximum(obs[idx])
Plots.plot(x_obs, lin_interp.(x_obs), color=:grey, ylabel="NMAD [m]", label="linear interpolation")
Plots.plot!(x_obs, exp_model(x_obs, fit.param), label="exponential fit")
Plots.plot!(x_obs, linear_model(x_obs, fit_lin.param), label="linear fit")
Plots.scatter!(bin_centers, nmad, xlabel="elevation [m]", ylabel="nmad", label="samples")

## choose a model for sigma
# σ_dh = exp_model(obs[idx], fit.param)
# σ_dh = linear_model(obs[idx], fit_lin.param)
σ_dh = lin_interp.(obs[idx])



### bias
# calculate the bias for each elevation bin
bias_bin = median.(dh_binned)

## model the bias as a function of elevation
lin_interp_bias = linear_interpolation(bin_centers, bias_bin, extrapolation_bc=Interpolations.Line())
@. bias_model(x,p) = -p[1]*exp(x*p[2])
fit_bias = LF.curve_fit(bias_model, bin_centers, bias_bin, [0.0,6e-4])
Plots.plot(x_obs, lin_interp_bias.(x_obs), color=:grey, ylabel="NMAD [m]", label="linear interpolation")
Plots.plot!(x_obs, bias_model(x_obs, fit_bias.param), label="exponential fit")
Plots.scatter!(bin_centers, bias_bin, xlabel="elevation [m]", ylabel="bias = mean difference", label="samples")

## chose a model for bias
# bias_dh  = bias_model(obs[idx], fit_bias.param)
bias_dh  = lin_interp_bias.(obs[idx])


### standardized, bias-corrected dh
z_dh = (dh  .- bias_dh ) ./ σ_dh

# plots after standardization
z_dh_binned, Z_dh_bin_centers = bin_equal_sample_size(obs[idx], z_dh, nsamples)
StatsPlots.violin(z_dh_binned, linewidth=0.5, color=:darkblue, legend=:false, fillalpha=0.6, ylabel="elevation difference [m]")
Plots.plot(
    StatsPlots.qqplot(StatsPlots.Normal(), dh, title="before standardization"),
    StatsPlots.qqplot(StatsPlots.Normal(), z_dh, title="after standardization \n and outlier filter")
)



# make table for GeoStats.jl
nx, ny = size(rec)
X = repeat(x, 1, ny)
X_data = reshape(X, nx*ny, 1)[idx]
Y = repeat(y, nx, 1)
Y_data = reshape(Y, nx*ny, 1)[idx]
table = (Z = z_dh, )
coords = [(first(X_data[i]), first(Y_data[i])) for i in 1:length(X_data)]

# compute empirical variogram
data = georef(table, coords)
U = data |> UniqueCoords()
gamma = EmpiricalVariogram(U, :Z; estimator=:matheron, nlags=700) #, maxlag=8e3)

Plots.scatter(gamma.abscissa, gamma.ordinate)

# fit the function for GaussianRandomFields manually
function make_custom_γ(Cor_fct...)
    n = length(Cor_fct)
    function custom_variogram(x, params)
        i_not_defined = findall(params .<= 0)
        params[i_not_defined] .= eps()
        # ep1 = Cor_fct[1](params[1], σ=params[2])
        # ep2 = Cor_fct[2](params[3], σ=params[4])
        # γ  = n .- GaussianRandomFields.apply.(ep1, x) .- GaussianRandomFields.apply.(ep2, x) # x = distance

        γ = 0
        for (i, f) in enumerate(Cor_fct)
            ep = f(params[(i-1)*2+1], σ=params[i*2])
            γ = γ .+ GaussianRandomFields.apply.(ep, x)
        end


        return γ
    end
    return custom_variogram
end

custom_model = make_custom_γ(GaussianRandomFields.Exponential, GaussianRandomFields.Exponential)

i_fit = findall(gamma.ordinate .!= 0) # .&& gamma.abscissa .< 1e6)           # some have values 0 at random distances, not sure why...
ff = LF.curve_fit(custom_model, gamma.abscissa[i_fit], 1 .- gamma.ordinate[i_fit], [1e4, 1.0, 4e6, 1.0])

Plots.plot([0; gamma.abscissa[i_fit]], custom_model([0; gamma.abscissa[i_fit]], ff.param))
Plots.scatter!(gamma.abscissa[i_fit], 1 .- gamma.ordinate[i_fit])


Plots.plot!([0; gamma.abscissa[i_fit]], GaussianRandomFields.apply.(Exponential(ff.param[1], σ=ff.param[2]), [0; gamma.abscissa[i_fit]]))
Plots.plot!([0; gamma.abscissa[i_fit]], GaussianRandomFields.apply.(Exponential(ff.param[3], σ=ff.param[4]), [0; gamma.abscissa[i_fit]]))


cov1 = CovarianceFunction(2, Exponential(ff.param[1], σ=ff.param[2]))
pts = range(0, step=900, length=101) # 1001 is the number of points
ptsy = range(0, step=900, length=1001) # 1001 is the number of points
grf = GaussianRandomField(cov1, Spectral(), pts, ptsy, n=5)
Plots.heatmap(grf)

# coarse_grid = round(ff.param[3], sigdigits=2)
coarse_grid = ff.param[1] / 5
tmpf = "coarse.nc"
coarse = gdalwarp(obs_aero_file; gr=coarse_grid, dest=tmpf)
idx_c = findall(vec(coarse .> 0.0))
nx_c, ny_c = size(coarse)
x_coarse = ncread(tmpf, "x")
X_c = repeat(x_coarse, 1, ny_c)
X_data_c = sort(reshape(X_c, nx_c*ny_c)[idx_c])
y_coarse = ncread(tmpf, "y")
Y_c = repeat(y_coarse, 1, nx_c)
Y_data_c = sort(reshape(Y_c, nx_c*ny_c)[idx_c])



cov = CovarianceFunction(2, Exponential(ff.param[3], σ=ff.param[4]))
grf = GaussianRandomField(cov, Spectral(), X_data_c, Y_data_c)
Plots.heatmap(grf)

