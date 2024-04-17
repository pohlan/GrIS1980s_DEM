using svd_IceSheetDEM
using NCDatasets, Interpolations, DataFrames, CSV, ProgressMeter, GeoStats, LsqFit, ParallelStencil, ImageMorphology, ImageSegmentation, JLD2
using StatsBase
import Plots
Plots.scalefontsizes(2.0)

# @init_parallel_stencil(Threads, Float64, 1)
# print("Number of threads:"); println(Threads.nthreads())

# stuff that's going to be in src or input to functions
gr = 1200
no_data_value = -9999.0
get_ix(i,nx) = i % nx == 0 ? nx : i % nx
get_iy(i,nx) = cld(i,nx)
get_global_i(ix, iy, nx) = nx * (iy-1) + ix
function replace_missing!(A,c)
    A[ismissing.(A)] .= c
    return A
end
# define paths
SEQ_output_path = "output/SEQ/"
fig_path        = joinpath(SEQ_output_path, "figures/")
sims_path       = joinpath(SEQ_output_path, "simulations/")
grimp_path      = "data/grimp/surface/"

mkpath(fig_path)
mkpath(sims_path)

# masks
ds_mask    = NCDataset("data/gris-imbie-1980/imbie_mask_g$(gr).nc")["Band1"][:,:]
bedm_mask  = NCDataset("data/bedmachine/bedmachine_g$(gr).nc")["mask"][:,:]
I_no_ocean = findall(.!ismissing.(vec(ds_mask)) .&& vec((bedm_mask .!= 1)))

# GrIMP (recent version)
ds_grimp = NCDataset(joinpath(grimp_path,"grimp_geoid_corrected_g$(gr).nc"))
grimp = ds_grimp["surface"][:,:]
x  = sort(ds_grimp["x"])
y  = sort(ds_grimp["y"])

# dhdt
dhdt = NCDataset("data/dhdt/CCI_GrIS_RA_SEC_5km_Vers3.0_2021-08-09_g$(gr)_1994-1996.nc")["Band1"][:,:]
replace_missing!(dhdt, 0)

# AERODEM
h_aero = NCDataset("data/aerodem/aerodem_g$(gr)_aligned.nc")["surface"][:,end:-1:1]
dh_aero_all = grimp - h_aero
idx_aero = findall(.!ismissing.(vec(dh_aero_all)) .&& .!ismissing.(vec(ds_mask)) .&& vec((bedm_mask .!= 1)))
df_aero  = DataFrame(:x       => x[get_ix.(idx_aero, length(x))],
                     :y       => y[get_iy.(idx_aero, length(x))],
                     :h_grimp => grimp[idx_aero],
                     :dh      => dh_aero_all[idx_aero],
                     :idx     => idx_aero,
                     :source .=> :aerodem )

# ATM
df_atm = CSV.read("data/grimp/surface/atm_new_grimp_interp_pts_no_registr.csv", DataFrame)
df_atm[!,:source] .= :atm
mindist = 1e4 # minimum distance between aerodem and atm points, ToDo: make this dependent on dhdt from the ATM instead
max_dhdt = std(dhdt[dhdt .!= 0])
function choose_atm(x_, y_)::Bool
    ix    = findmin(abs.(x_ .- x))[2]
    iy    = findmin(abs.(y_ .- y))[2]
    iglob = get_global_i(ix, iy, length(x))
    # filter values outside the ice sheet
    is_in_icesheet = iglob ∈ I_no_ocean
    # filter values overlapping or close to aerodem
    dist_to_aero = minimum(pairwise(Distances.euclidean, [x_ y_], [df_aero.x df_aero.y], dims=1)[:])
    is_above_mindist = dist_to_aero > mindist
    # filter values with high absolute dhdt
    has_low_dhdt = 0.0 < abs(dhdt[ix,iy]) < max_dhdt
    return is_in_icesheet & is_above_mindist & has_low_dhdt
end
println("selecting flightline values...")
id = sort(StatsBase.sample(1:size(df_atm,1), 100000, replace=false))  # without pre-selecting a subsample of points the filtering takes a long time
keepat!(df_atm, id)
filter!([:x,:y] => choose_atm, df_atm)
# already remove some outliers here, improves standardization
atm_to_delete = findall(abs.(df_atm.dh) .> 5 .* mad(df_atm.dh))
deleteat!(df_atm, atm_to_delete)

# merge aerodem and atm data
df_all = vcat(df_aero, df_atm, cols=:intersect)

# plot
dh_plot = Array{Union{Float32,Missing}}(missing, (length(x),length(y)))
dh_plot[df_aero.idx] = df_aero.dh
Plots.heatmap(x,y, dh_plot', cmap=:RdBu, clims=(-20,20))
Plots.scatter!(df_atm.x, df_atm.y, marker_z=df_atm.dh, color=:RdBu, markersize=0.5, markerstrokewidth=0, legend=false)
Plots.savefig(joinpath(fig_path,"data_non-standardized.png"))

# mean trend as a fct of elevation
dh_binned, ELV_bin_centers = svd_IceSheetDEM.bin_equal_sample_size(df_all.h_grimp, df_all.dh, 12000)  # 16000
itp_bias = linear_interpolation(ELV_bin_centers, mean.(dh_binned),  extrapolation_bc = Interpolations.Flat())
Plots.scatter(ELV_bin_centers, mean.(dh_binned))
ELV_plot = minimum(df_all.h_grimp):10:maximum(df_all.h_grimp)
Plots.plot!(ELV_plot, itp_bias(ELV_plot), label="linear interpolation", legend=:bottomright)
Plots.savefig(joinpath(fig_path,"bias_vs_elevation.png"))

# standard deviation as a fct of elevation
itp_std = linear_interpolation(ELV_bin_centers, mad.(dh_binned),  extrapolation_bc = Interpolations.Flat())
Plots.scatter(ELV_bin_centers, mad.(dh_binned))
ELV_plot = minimum(df_all.h_grimp):10:maximum(df_all.h_grimp)
Plots.plot!(ELV_plot, itp_std(ELV_plot), label="linear interpolation", legend=:topright)
Plots.savefig(joinpath(fig_path,"mad_vs_elevation.png"))

# standardize
df_all.dh_detrend   = (df_all.dh .- itp_bias(df_all.h_grimp)) ./ itp_std(df_all.h_grimp)
# remove outliers after standardizing
all_to_delete = findall(abs.(df_all.dh_detrend) .> 7 .* mad(df_all.dh_detrend))
deleteat!(df_all, all_to_delete)
# divide again by std as dividing by mad doesn't normalize fully
std_dh_detrend      = std(df_all.dh_detrend)
df_all.dh_detrend ./= std_dh_detrend
# plot standardized histograms
Plots.histogram(df_all.dh_detrend, label="Standardized observations", xlims=(-10,10), normalize=:pdf, nbins=1000, wsize=(600,500), linecolor=nothing)
Plots.histogram!((df_all.dh .- mean(df_all.dh)) ./ std(df_all.dh), label="Standardized without binning", normalize=:pdf, linecolor=nothing)
Plots.plot!(Normal(), lw=1, label="Normal distribution", color="black")
Plots.savefig(joinpath(fig_path,"histogram_standardization.png"))

# plot again after standardizing
Plots.scatter(df_all.x, df_all.y, marker_z=df_all.dh_detrend, label="", markersize=2.0, markerstrokewidth=0, cmap=:RdBu, clims=(-4,4), aspect_ratio=1, xlims=(-7e5,8e5), xlabel="Easting [m]", ylabel="Northing [m]", colorbar_title="[m]", title="Standardized elevation difference (GrIMP - historic)", grid=false, wsize=(1700,1800))
# Plots.scatter(df_varg.x, df_varg.y, marker_z=df_varg.dh_detrend, label="", markersize=2.0, markerstrokewidth=0, cmap=:RdBu, clims=(-2,2), aspect_ratio=1, xlims=(-7e5,8e5), xlabel="Easting [m]", ylabel="Northing [m]", colorbar_title="[m]", title="Standardized elevation difference (GrIMP - historic)", grid=false, wsize=(1700,1800))
Plots.savefig(joinpath(fig_path,"data_standardized.png"))

# variogram
println("Variogram...")
table_all  = (; Z=df_all.dh_detrend)
coords_all = [(df_all.x[i], df_all.y[i]) for i in 1:length(df_all.x)]
data       = georef(table_all,coords_all)
gamma      = EmpiricalVariogram(data, :Z; nlags=400,  maxlag=7e5, estimator=:cressie)
# fit a covariance function
function custom_var(x, params)
    if any(params[[1,3,5]] .>= 8e5)
        return -9999.0 .* ones(length(x))
    end
    γ1 = SphericalVariogram(range=params[1], sill=params[2])      # zero nugget
    γ2 = SphericalVariogram(range=params[3], sill=params[4])
    γ3 = SphericalVariogram(range=params[5], sill=params[6])
    f = γ1.(x) .+ γ2.(x) .+ γ3.(x)
    return f
end
ff = LsqFit.curve_fit(custom_var, gamma.abscissa, gamma.ordinate, [1e5, 0.2, 1e4, 0.2, 5e5, 0.5]);
varg = SphericalVariogram(range=ff.param[1], sill=ff.param[2]) +
       SphericalVariogram(range=ff.param[3], sill=ff.param[4]) +
       SphericalVariogram(range=ff.param[5], sill=ff.param[6])
Plots.scatter(gamma.abscissa, gamma.ordinate, label="all observations", xscale=:log10)
Plots.plot!(gamma.abscissa,custom_var(gamma.abscissa, ff.param), label="LsqFit fit", lw=2)
Plots.savefig(joinpath(fig_path,"variogram.png"))

# do simulations
maxn   = 50     # maximum neighbours taken into account
n_sims = 20     # number of simulations

table_input  = (; Z=Float32.(df_all.dh_detrend))
coords_input = [(xi,yi) for (xi,yi) in zip(df_all.x,  df_all.y)]
geotable_input   = georef(table_input,coords_input)

# output as pointset rather than grid because grid includes a lot of unnecessary points in the ocean etc, makes it a lot slower
I_no_ocean_no_aero = findall(.!ismissing.(vec(ds_mask)) .&& vec((bedm_mask .!= 1)) .&& ismissing.(vec(dh_aero_all)))
xnew = x[get_ix.(I_no_ocean_no_aero, length(x))]
ynew = y[get_iy.(I_no_ocean_no_aero, length(x))]
grid_output = PointSet([Point(xi,yi) for (xi,yi) in zip(xnew, ynew)])

process = GaussianProcess(varg)
method = SEQMethod(maxneighbors=maxn)

println("Do SEQ simulations...")
tic = Base.time()
sims = rand(process, grid_output, geotable_input, n_sims, method)
toc = Base.time() - tic
print("Time for simulations in minutes: "); println(toc/60)

# # for point mesh
# Plots.scatter(xnew, ynew, marker_z=sim.Z, markerstrokewidth=0, markersize=1.0, cmap=:RdBu, clims=(-2,2))
# Plots.savefig("interpol.png")

# @parallel_indices ind function do_seq!(ms)
#     if 1 <= ind <= size(ms,2)
#         sim = rand(process, grid_output, geotable_input, method)
#         @views ms[:,ind] = sim.Z
#     end
#     return
# end
# ms = zeros(length(grid_output), n_sims)
# tic = Base.time()
# @parallel do_seq!(ms)
# toc = Base.time() - tic

## destandardize and move to higher resolution
grimp_tg = copy(grimp)
grimp_tg[ismissing.(grimp)] .= 0

# read in data at higher res
gr_up                   = 600
h_aero_g600             = NCDataset("data/aerodem/aerodem_g$(gr_up)_aligned.nc")["surface"][:,end:-1:1]
ds_mask_g600            = NCDataset("data/gris-imbie-1980/imbie_mask_g$(gr_up).nc")["Band1"][:,:]
bedm_mask_g600          = NCDataset("data/bedmachine/bedmachine_g$(gr_up).nc")["mask"][:,:]
I_no_ocean_no_aero_g600 = findall(.!ismissing.(vec(ds_mask_g600)) .&& (vec(bedm_mask_g600) .!= 1) .&& ismissing.(vec(h_aero_g600)))

for (i,s) in enumerate(sims)
    dest_file_0  = joinpath(sims_path,"SEQ_maxngh_$(maxn)_g$(gr)_id_$(i).nc")
    dest_file_up = joinpath(sims_path,"SEQ_maxngh_$(maxn)_g$(gr_up)_id_$(i)_upsampled.nc")

    h_predict_all                        = zeros(size(h_aero))
    h_predict_all[.!ismissing.(h_aero)] .= h_aero[.!ismissing.(h_aero)]
    h_predict_all[I_no_ocean_no_aero]    = grimp_tg[I_no_ocean_no_aero] .- (s.Z .*std_dh_detrend.*itp_std.(grimp_tg[I_no_ocean_no_aero]) .+ itp_bias.(grimp_tg[I_no_ocean_no_aero]))
    h_predict_all[h_predict_all .<= 0 .|| isnan.(h_predict_all)] .= no_data_value
    svd_IceSheetDEM.save_netcdf(dest_file_0, joinpath(grimp_path,"grimp_geoid_corrected_g$(gr).nc"), [h_predict_all], ["surface"], Dict("surface" => Dict{String, Any}()))

    # moving to higher resolution (upsampling)
    upsampled_temp = joinpath(sims_path,"SEQ_upsampled_temp.nc")
    gdalwarp(dest_file_0; gr=gr_up, srcnodata="-9999", dstnodata="-9999", dest=upsampled_temp)
    h_predict_warped = NCDataset(upsampled_temp)["Band1"][:,:]
    rm(upsampled_temp)
    h_predict_warped[ismissing.(h_predict_warped)] .= 0
    # put together in new array
    h_new = zeros(size(h_aero_g600))
    h_new[.!ismissing.(h_aero_g600)] .= h_aero_g600[.!ismissing.(h_aero_g600)]     # use original observations where available
    h_new[I_no_ocean_no_aero_g600]   .= h_predict_warped[I_no_ocean_no_aero_g600]  # use upsampled SEQ simulation to fill the interior
    h_new[h_new .<= 0 .|| isnan.(h_new)] .= no_data_value
    svd_IceSheetDEM.save_netcdf(dest_file_up, joinpath(grimp_path,"grimp_geoid_corrected_g$(gr_up).nc"), [h_new], ["surface"], Dict("surface" => Dict{String, Any}()))
end
