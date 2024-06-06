using svd_IceSheetDEM
using NCDatasets, Interpolations, DataFrames, CSV, ProgressMeter, GeoStats, JLD2
using StatsBase, Distributions
import Plots
Plots.scalefontsizes(2.0)

# stuff that's going to be in src or input to functions
gr = 1200
no_data_value = -9999.0
# get_ix(i,nx) = i % nx == 0 ? nx : i % nx    # use LinearIndices instead
# get_iy(i,nx) = cld(i,nx)
# get_global_i(ix, iy, nx) = nx * (iy-1) + ix

function nanfun(f,a)
    id_numbers = findall(.!isnan.(a))
    return f(a[id_numbers])
end
# define paths
SEQ_output_path = "output/SEQ/"
fig_path        = joinpath(SEQ_output_path, "figures/")
sims_path       = joinpath(SEQ_output_path, "simulations/")

mkpath(fig_path)
mkpath(sims_path)

# ToDo: get filenames from functions as in main.jl
rec_file = "output/rec_files/rec_lambda_1e7_g600_r200.nc"
bedm_file = "data/bedmachine/bedmachine_g600.nc"
obs_aero_file = "data/aerodem/aerodem_rm-filtered_geoid-corr_g600.nc"
obs_ATM_file = get_atm_file()
dhdt_file = "data/dhdt/CCI_GrIS_RA_SEC_5km_Vers3.0_2021-08-09_g$(gr)_1994-1996.nc"
mask_file = "data/gris-imbie-1980/imbie_mask_g$(gr).nc"

dh_atm_file = joinpath(dirname(obs_ATM_file), "GrIMPv1_minus_atm.csv")

# define variogram function to fit
custom_var(params) = SphericalVariogram(range=params[1], sill=params[4] ./ sum(params[4:6]), nugget=params[7]) +
                     SphericalVariogram(range=params[2], sill=params[5] ./ sum(params[4:6]), nugget=params[8]) +
                     SphericalVariogram(range=params[3], sill=params[6] ./ sum(params[4:6]), nugget=params[9])
param_cond(params) = all(params[1:3] .< 1e6) && all(params[7:9] .> 0)
# initial guess for parameters
p0 = [1.2e4, 1.3e5, 6e5, 0.1, 0.3, 0.6, 0.,0.,0.]
_, varg, _, destand, I_no_ocean, _ = prepare_random_sims(rec_file, bedm_file, obs_aero_file, obs_ATM_file, dhdt_file, mask_file;
                                                       dh_atm_file, fig_path, custom_var, param_cond, p0, nbins1=6, nbins2=10)

# do simulations
n_sims = 10     # number of simulations

table_input    = (; Z=Float32.(df_all.dh_detrend))
coords_input   = [(xi,yi) for (xi,yi) in zip(df_all.x,  df_all.y)]
geotable_input = georef(table_input,coords_input)

# subsample for longer ranges to have an effect (limited by maxneighbors parameter in SEQ)
frac_subs = 4
# id_all = Vector(1:size(df_all,1))
# id = 1:(length(id_all)-cld(length(id_all),frac_subs))
# df_SEQ = df_all[id,:]
# df_test = df_all[setdiff(id_all,id),:]
# filter!( :source => x -> x==:atm, df_test)
Plots.scatter(df_SEQ.x, df_SEQ.y, marker_z=df_SEQ.dh_detrend, label="", markersize=2.0, markerstrokewidth=0, cmap=:RdBu, clims=(-4,4), aspect_ratio=1, xlims=(-7e5,8e5), xlabel="Easting [m]", ylabel="Northing [m]", colorbar_title="[m]", title="Standardized elevation difference (GrIMP - historic)", grid=false, wsize=(1700,1800))
Plots.scatter(df_test.x, df_test.y, marker_z=df_test.dh_detrend, label="", markersize=2.0, markerstrokewidth=0, cmap=:RdBu, clims=(-4,4), aspect_ratio=1, xlims=(-7e5,8e5), xlabel="Easting [m]", ylabel="Northing [m]", colorbar_title="[m]", title="Standardized elevation difference (GrIMP - historic)", grid=false, wsize=(1700,1800))
Plots.savefig(joinpath(fig_path,"data_sampled_for_SEQ.png"))


# simulate data on higher resolution
gr_up                   = 1200
h_aero_g600             = NCDataset("data/aerodem/aerodem_g$(gr_up)_aligned.nc")["surface"][:,:]
h_aero_g600             = h_aero_g600[:,end:-1:1]
ds_mask_g600            = NCDataset("data/gris-imbie-1980/imbie_mask_g$(gr_up).nc")["Band1"][:,:]
bedm_mask_g600          = NCDataset("data/bedmachine/bedmachine_g$(gr_up).nc")["mask"][:,:]
I_no_ocean_no_aero_g600 = findall(.!ismissing.(vec(ds_mask_g600)) .&& (vec(bedm_mask_g600) .!= 1) .&& ismissing.(vec(h_aero_g600)))
x_g600 = NCDataset("data/aerodem/aerodem_g$(gr_up)_aligned.nc")["x"][:]
y_g600 = NCDataset("data/aerodem/aerodem_g$(gr_up)_aligned.nc")["y"][:]

# output as pointset rather than grid because grid includes a lot of unnecessary points in the ocean etc, makes it a lot slower
I_no_ocean_no_aero = findall(.!ismissing.(vec(ds_mask)) .&& vec((bedm_mask .!= 1)) .&& ismissing.(vec(dh_aero_all)))
xnew = x_g600[get_ix.(I_no_ocean_no_aero_g600, length(x_g600))]
ynew = y_g600[get_iy.(I_no_ocean_no_aero_g600, length(x_g600))]
n1   = 1*10^4
nend = 2*10^5
# grid_output = PointSet([Point(xi,yi) for (xi,yi) in zip(xnew[n1:nend], ynew[n1:nend])])
# grid_output = PointSet([Point(xi,yi) for (xi,yi) in zip(xnew, ynew)])
grid_output = PointSet([Point(xi,yi) for (xi,yi) in zip(x[get_ix.(I_no_ocean, length(x))], y[get_iy.(I_no_ocean, length(x))])])

process = GaussianProcess(varg)

maxn   = 1000    # maximum neighbours taken into account
method = SEQMethod(maxneighbors=maxn)


# # aerodem geotable
# id_aero = findall(df_all.source .== :aerodem)
# geotable_aero  = view(geotable_input, id_aero)
# # atm geotable
# id_atm = findall(df_all.source .== :atm)
# geotable_atm   = view(geotable_input, id_atm)
# df_aero.dh_detrend   = (df_aero.dh .- itp_bias(df_aero.h_grimp)) ./ itp_std(df_aero.h_grimp)
# df_aero.dh_detrend ./= std_dh_detrend
# df_sgs = vcat(df_all, df_aero, cols=:intersect)
# table_sgs    = (; Z=Float32.(df_sgs.dh_detrend))
# coords_sgs   = [(xi,yi) for (xi,yi) in zip(df_sgs.x,  df_sgs.y)]
# geotable_input = georef(table_sgs,coords_sgs)




# partitioning in Meshes.jl: https://juliageometry.github.io/MeshesDocs/stable/algorithms/partitioning.html
# prts = partition(PointSet(Point.(coords_input)), BlockPartition(4e5,neighbors=true))
# Meshes.metadata(prts)[:neighbors]
flds = folds(geotable_input, BlockFolding(4e5))
# inspired by https://github.com/JuliaEarth/GeoStatsValidation.jl/blob/main/src/cverrors/wcv.jl#L75 but there not implemented for random simulations conditioned on data
dif_blocks     = [Float64[] for i in flds]
NNdist_blocks  = [Float64[] for i in flds]
xcoord_blocks  = [Float64[] for i in flds]
ycoord_blocks  = [Float64[] for i in flds]
p = Plots.plot()
# fig, ax, hm = viz(geotable_input.geometry[1])
@showprogress for (j,fs) in enumerate(flds)
    # find the neighbors that the folds routine (https://github.com/JuliaEarth/GeoStatsBase.jl/blob/master/src/folding/block.jl) leaves out
    # there might be a mistake in the partitioning routine in Meshes.jl, the neighbors don't make sense (also not tested well)
    neighbors = Vector(1:length(geotable_input.Z))
    deleteat!(neighbors, unique(sort([fs[1];fs[2]])))
    append!(fs[1],neighbors)

    # sdat  = vcat(view(geotable_input, fs[1]), geotable_aero)
    sdat  = view(geotable_input, fs[1])
    stest = view(domain(geotable_input), fs[2])
    @assert length(sdat.Z) > length(stest)

    # subsample again
    # id_subspl = sort(StatsBase.sample(1:length(sdat.Z), ceil(Int,length(sdat.Z)*frac_subs); replace=false))
    # id = 1:(length(sdat.Z)-cld(length(sdat.Z),frac_subs))
    # sdat = view(sdat, id_subspl)
    # if j==8

    # simulations
    sims_test    = rand(process, stest, sdat, 20, method)
    var_sims     = zeros(length(sims_test))
    test_data    = view(geotable_input, fs[2]).Z
    dif_blocks[j] = mean(sims_test).Z .- test_data


    # save distances
    # ids_nn  = zeros(Int,length(stest))
    # dist_nn = zeros(length(stest))
    # for (i,st) in enumerate(stest)
    #     id_st, dist = searchdists(st, KNearestSearch(domain(sdat),1))
    #     ids_nn[i]  = id_st[1]
    #     dist_nn[i] = dist[1]
    # end
    # @assert length(dist_nn) == length(test_data)
    # NNdist_blocks[j] = dist_nn

    # # # save coordinates
    crds = coordinates.(stest)
    xcoord_blocks[j] = first.(crds)
    ycoord_blocks[j] = last.(crds)

    # plot
        # Plots.scatter!(df_all.x[fs[1]],df_all.y[fs[1]], color="grey",markersize=2, markerstrokewidth=0)
        # Plots.scatter!(df_all.x[fs[2]],df_all.y[fs[2]], color=cols[j], markerstrokewidth=0,label=string(j))
        # subd = view(sdat, sort(unique(ids_nn)))
        # b=[m.coords for m in subd.geometry]
        # Plots.scatter!(first.(b), last.(b),     color="black", markerstrokewidth=0, label=string(j))
        # Plots.scatter!(first.(crds),last.(crds), marker_z=mean(sims_test).Z, cmap=:bwr, clims=(-3,3), markerstrokewidth=0)
        # Plots.scatter!(first.(crds),last.(crds), marker_z=test_data, cmap=:bwr, clims=(-3,3), markerstrokewidth=0)

        # viz!(domain(sdat), color=cols[j])
        # viz!(domain(subd), color="black", markersize=10)
        # viz!(stest, color=mean(sims_test).Z)
    # if j == 6
    #     break
    # end
end
Plots.plot(p)

dists = vcat(NNdist_blocks...)
difs  = vcat(dif_blocks...)
xcrds = vcat(xcoord_blocks...)
ycrds = vcat(ycoord_blocks...)
# Plots.scatter(dists,difs,markerstrokewidth=0,markersize=1,color="grey")
dif_binned, dist_bin_centers = svd_IceSheetDEM.bin_equal_sample_size(dists[1:length(difs)], difs, 6000)  # 16000
Plots.scatter(dist_bin_centers./1e3, nanfun.(median,dif_binned), title="maxn=$maxn", xlabel="distance to closest neighbor", ylabel="median", label="", xscale=:log10, ylims=(-0.3,0.2))
Plots.scatter(dist_bin_centers./1e3, mad.(dif_binned), title="maxn=$maxn", xlabel="distance to closest neighbor", ylabel="nmad", label="", xscale=:log10)
Plots.savefig("median.png")

Plots.scatter(xcrds,ycrds,marker_z=difs, cmap=:bwr, clims=(-5,5),markerstrokewidth=0,markersize=2)
Plots.scatter(xcrds, ycrds, marker_z=geotable_input.Z, cmap=:bwr, clims=(-5,5),markerstrokewidth=0,markersize=2)

nm = 5
Plots.scatter(vcat(setdiff(xcoord_blocks,nm)...),vcat(setdiff(ycoord_blocks,nm)...), markerstrokewidth=0)
Plots.scatter!(xcoord_blocks[nm],ycoord_blocks[nm], markerstrokewidth=0)
Plots.scatter!(xcoord_blocks[nm],ycoord_blocks[nm], markerstrokewidth=0)


# println("Do SEQ simulations...")
# tic = Base.time()
# Random.seed!(1234)
# sims = rand(process, grid_output, geotable_input, n_sims, method)
# toc = Base.time() - tic
# print("Time for simulations in minutes: "); println(toc/60)


# test
test_output = PointSet([Point(xi,yi) for (xi,yi) in zip(df_test.x, df_test.y)])
sims_test = rand(process, test_output, geotable_input, 30, method)
err = []
dd  = zeros(length(sims_test), length(sims_test[1].Z))
for (ism,sm) in enumerate(sims_test)
    push!(err,median(abs.(sm.Z .- df_test.dh_detrend)))
    dd[ism,:] = sm.Z .- df_test.dh_detrend
end
println(median(err))
m_per_p = median(dd,dims=1)[:]
# st_per_p = std(dd,dims=1)[:]
Plots.scatter(df_test.x,df_test.y,marker_z=m_per_p, cmap=:bwr, clims=(-2,2), markerstrokewidth=0, markersize=1, aspect_ratio=1, title="maxn=$maxn", label="")
# Plots.scatter(df_test.x,df_test.y,marker_z=st_per_p, markerstrokewidth=0, markersize=1, aspect_ratio=1, title="maxn=$maxn",label="")


xn1, xnn = extrema(xnew[n1:nend])
yn1, ynn = extrema(ynew[n1:nend])
function choose_bla(x_, y_)::Bool
    ix    = xn1 < x_ < xnn
    iy    = yn1 < y_ < ynn
    return ix && iy
end
df_plot = filter([:x,:y] => choose_bla, df_SEQ)
Plots.scatter(df_plot.x, df_plot.y, marker_z=df_plot.dh_detrend, markerstrokewidth=0, markersize=4.0, cmap=:RdBu, clims=(-2.5,2.5), aspect_ratio=1, xlabel="Easting [m]", ylabel="Northing [m]", colorbar_title="[m]", grid=false, wsize=(900,1000))
Plots.scatter(xnew[n1:nend], ynew[n1:nend], marker_z=sims[1].Z, markerstrokewidth=0, markersize=1.0, cmap=:RdBu,  clims=(-2.5,2.5), aspect_ratio=1, xlabel="Easting [m]", ylabel="Northing [m]", colorbar_title="[m]", grid=false, wsize=(900,1000), title="maxn=$maxn, frac_subs=$frac_subs")
# Plots.scatter(xnew[n1:nend], ynew[n1:nend], marker_z=mean(sims).Z, markerstrokewidth=0, markersize=1.0, cmap=:RdBu,  clims=(-2.5,2.5), aspect_ratio=1, xlabel="Easting [m]", ylabel="Northing [m]", colorbar_title="[m]", grid=false, wsize=(900,1000), title="maxn=$maxn, frac_subs=$frac_subs")
# Plots.savefig(joinpath(fig_path,"interpolation_simulation_maxn_$(maxn)_frac_$(frac_subs)_aeroonly.png"))
Plots.scatter(x[get_ix.(I_no_ocean, length(x))], y[get_iy.(I_no_ocean, length(x))], marker_z=sims[1].Z, markerstrokewidth=0, markersize=1.0, cmap=:RdBu,  clims=(-2.5,2.5), aspect_ratio=1, xlabel="Easting [m]", ylabel="Northing [m]", colorbar_title="[m]", grid=false, wsize=(900,1000), title="maxn=$maxn, frac_subs=$frac_subs")
Plots.scatter(df_SEQ.x, df_SEQ.y, marker_z=df_SEQ.dh_detrend, markerstrokewidth=0, markersize=4.0, cmap=:RdBu, clims=(-2.5,2.5), aspect_ratio=1, xlabel="Easting [m]", ylabel="Northing [m]", colorbar_title="[m]", grid=false, wsize=(900,1000))


table_     = (; Z=sims[1].Z)
# coords_    = [(xnew[i], ynew[i]) for i in n1:nend]
coords_    = [(x[get_ix.(i, length(x))], y[get_iy.(i, length(x))]) for i in I_no_ocean]
data       = georef(table_,coords_)
gamma      = EmpiricalVariogram(data, :Z; nlags=200,  maxlag=1e5, estimator=:cressie)
ff = LsqFit.curve_fit(custom_var, gamma.abscissa, gamma.ordinate, [1e5, 0.2, 1e4, 0.2, 5e5, 0.5]);


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
grimp_tg = NCDataset("data/bedmachine/bedmachine_g$(gr_up).nc")["surface"][:,:]
grimp_tg[ismissing.(grimp_tg)] .= 0


function destandardize_and_save(zvals, dest_file_0, dest_file_up)
    h_predict_all                        = zeros(size(h_aero_g600))
    h_predict_all[.!ismissing.(h_aero_g600)] .= h_aero_g600[.!ismissing.(h_aero_g600)]
    h_predict_all[I_no_ocean_no_aero_g600]    = grimp_tg[I_no_ocean_no_aero_g600] .- (zvals .*std_dh_detrend.*itp_std.(grimp_tg[I_no_ocean_no_aero_g600]) .+ itp_bias.(grimp_tg[I_no_ocean_no_aero_g600]))
    h_predict_all[h_predict_all .<= 0 .|| isnan.(h_predict_all)] .= no_data_value
    svd_IceSheetDEM.save_netcdf(dest_file_0, joinpath(grimp_path,"grimp_geoid_corrected_g$(gr_up).nc"), [h_predict_all], ["surface"], Dict("surface" => Dict{String, Any}()))

    dh_aero_all_g600 = grimp_tg - h_aero_g600
    dh_to_plot = zeros(size(dh_aero_all_g600))
    ifll = findall(.!ismissing.(vec(dh_aero_all_g600)) .&& .!ismissing.(vec(ds_mask_g600)) .&& vec((bedm_mask_g600 .!= 1)))
    dh_to_plot[ifll] .= (dh_aero_all_g600[ifll] .- itp_bias(grimp_tg[ifll])) ./ itp_std(grimp_tg[ifll])
    dh_to_plot ./= std_dh_detrend
    dh_to_plot[I_no_ocean_no_aero_g600] .= zvals
    Plots.heatmap(x_g600, sort(y_g600), dh_to_plot', label="", cmap=:RdBu, clims=(-4,4), aspect_ratio=1, xlims=(-7e5,8e5), xlabel="Easting [m]", ylabel="Northing [m]", colorbar_title="[m]", title="Standardized elevation difference (GrIMP - historic)", grid=false, wsize=(1700,1800))
    Plots.savefig(joinpath(fig_path,"interpolated_standardized.png"))

    # moving to higher resolution (upsampling)
    # upsampled_temp = joinpath(sims_path,"SEQ_upsampled_temp.nc")
    # gdalwarp(dest_file_0; gr=gr_up, srcnodata="-9999", dstnodata="-9999", dest=upsampled_temp)
    # h_predict_warped = NCDataset(upsampled_temp)["Band1"][:,:]
    # rm(upsampled_temp)
    # h_predict_warped[ismissing.(h_predict_warped)] .= 0
    # # put together in new array
    # h_new = zeros(size(h_aero_g600))
    # h_new[.!ismissing.(h_aero_g600)] .= h_aero_g600[.!ismissing.(h_aero_g600)]     # use original observations where available
    # h_new[I_no_ocean_no_aero_g600]   .= h_predict_warped[I_no_ocean_no_aero_g600]  # use upsampled SEQ simulation to fill the interior
    # h_new[h_new .<= 0 .|| isnan.(h_new)] .= no_data_value
    # svd_IceSheetDEM.save_netcdf(dest_file_up, joinpath(grimp_path,"grimp_geoid_corrected_g$(gr_up).nc"), [h_new], ["surface"], Dict("surface" => Dict{String, Any}()))
    return
end

for (i,s) in enumerate(sims)
    zvals = s.Z
    dest_file_0  = joinpath(sims_path,"SEQ_maxngh_$(maxn)_g$(gr_up)_id_$(i).nc")
    dest_file_up = joinpath(sims_path,"SEQ_maxngh_$(maxn)_g$(gr_up)_id_$(i)_upsampled.nc")
    destandardize_and_save(zvals, dest_file_0, dest_file_up)
end

# save mean (should converge towards simple kriging solution)
zvals = mean(sims).Z
dest_file_0  = joinpath(sims_path,"SEQ_maxngh_$(maxn)_g$(gr)_mean.nc")
dest_file_up = joinpath(sims_path,"SEQ_maxngh_$(maxn)_g$(gr_up)_mean_upsampled.nc")
destandardize_and_save(zvals, dest_file_0, dest_file_up)
