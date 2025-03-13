using svd_IceSheetDEM, NetCDF, NCDatasets, Meshes, Statistics, StatsBase, GeoStats, CSV, DataFrames, Printf, JLD2, Distributions, ImageMorphology, StatsPlots, LaTeXStrings, Interpolations, UnPack
import Plots

# set target directories
main_output_dir = joinpath("output","validation")
mkpath(main_output_dir)
fig_dir         = joinpath(main_output_dir, "figures")
mkpath(fig_dir)

# for running the script interactively
# ARGS = [
#         "--shp_file", "data/gris-imbie-1980/gris-outline-imbie-1980_updated.shp",
#         "--grid_size", 600.0]

parsed_args         = parse_commandline(ARGS)
outline_shp_file    = parsed_args["shp_file"]
grd                 = parsed_args["grid_size"]
maxn                = parsed_args["maxn"]

################################################
# Preprocessing, standardization and variogram #
################################################

csv_preprocessing, jld2_preprocessing = svd_IceSheetDEM.prepare_obs(grd, outline_shp_file)

# get I_no_ocean, (de-)standardization functions and variogram from pre-processing
df_all = CSV.read(csv_preprocessing, DataFrame)
dict   = load(jld2_preprocessing)
@unpack I_no_ocean, idx_aero, gamma, href_file = dict
@unpack destandardize = svd_IceSheetDEM.get_stddization_fcts(jld2_preprocessing)
# varg = svd_IceSheetDEM.custom_var(params)
varg = svd_IceSheetDEM.get_var(gamma)

# make geotable
geotable = svd_IceSheetDEM.make_geotable(df_all.dh_detrend, df_all.x, df_all.y)

##############################################
# Standard cross-validation, leave block out #
##############################################

# create sets of training and test data
ℓ    = 2e5
flds = folds(geotable, BlockFolding(ℓ))

function evaluate_fun(i_train,i_test)
    interp = svd_IceSheetDEM.do_kriging(view(domain(geotable),i_test), view(geotable,i_train), varg; maxn)
    return interp.Z
end
println("Kriging cross-validation...")
difs, xc, yc = svd_IceSheetDEM.step_through_folds(flds, evaluate_fun, geotable, save_coords=true, save_distances=false)

# get indices
idxs = [Int[] for i in flds]
for (i,fs) in enumerate(flds)
    idxs[i] = fs[2]
end
idx = vcat(idxs...)

logℓ    = round(log10(ℓ),digits=1)
dest    = get_cv_file_kriging(grd, logℓ, maxn)
cv_dict = (; maxn, difs, xc, yc, idx, binfield1=df_all.bfield_1[idx], h_ref=df_all.h_ref[idx], grd, method="kriging")
jldsave(dest; cv_dict...)

##############################################

# plot map of difs
p = Plots.scatter(xc,yc,marker_z=difs, cmap=:bwr, clims=(-4,4), markersize=0.7, markerstrokewidth=0, label="", aspect_ratio=1, size=(500,700))
for (j,fs) in enumerate(flds)
    stest = view(domain(geotable), fs[2])
    ch    = GeoStats.convexhull(stest).rings[1].vertices.data
    xx    = [ustrip.(a.coords.x) for a in ch]
    yy    = [ustrip.(a.coords.y) for a in ch]
    Plots.plot!(p, xx,yy, lw=1, color="black", label="")
end
Plots.plot(p, size=(500,700))
Plots.savefig(joinpath(fig_dir, "map_validation_kriging_maxn$maxn.png"))

# uncertainty estimation
dem_ref                = NCDataset(href_file)["surface"][:,:]
dif_destd              = destandardize(cv_dict.difs, cv_dict.binfield1, cv_dict.h_ref, add_mean=false)
dh_binned, bin_centers = svd_IceSheetDEM.bin_equal_bin_size(cv_dict.h_ref, dif_destd, 14)
sitp, std_uncertainty  = uncertainty_from_cv(dh_binned, bin_centers, dem_ref)
dest_file              = get_std_uncrt_file(cv_dict.method, grd)
svd_IceSheetDEM.save_netcdf(dest_file, href_file, [std_uncertainty], ["std_uncertainty"], Dict("std_uncertainty" => Dict{String,Any}()))

#######################################################
# Determine good maxns through morphological gradient #
#######################################################
bedmachine_original, bedm_file    = svd_IceSheetDEM.create_bedmachine_grid(grd)
reference_file_g150, _, ref_file  = svd_IceSheetDEM.create_grimpv2(grd, dict["coreg_grid"], bedmachine_original)

# derive indices for cells to interpolate, this time interpolate the points we don't know, technically no validation!
ir_sim      = setdiff(I_no_ocean, idx_aero)  # indices that are in I_no_ocean but not in idx_aero
xsp = 650:1300
ysp = 250:800
x           = NCDataset(ref_file)["x"][:]
y           = NCDataset(ref_file)["y"][:]
ix = get_ix.(ir_sim, length(x))
iy = get_iy.(ir_sim, length(x))
ir_keep = findall(xsp[1].<ix.<xsp[end] .&& ysp[1].<iy.<ysp[end])
grid_output = PointSet([Point(xi,yi) for (xi,yi) in zip(x[ix[ir_keep]], y[iy[ir_keep]])])

ir_keep_df = findall(x[xsp[1]].<df_all.x.<x[xsp[end]] .&& y[ysp[1]].<df_all.y.<y[ysp[end]])
geotable = svd_IceSheetDEM.make_geotable(df_all.dh_detrend[ir_keep_df], df_all.x[ir_keep_df], df_all.y[ir_keep_df])

# loop through maxns and calculate the morphological gradient each time
maxns = [10, 100, 1500, 3000]
∑grad = zeros(length(maxns))
wallt = zeros(length(maxns))
m_interps = zeros(length(xsp), length(ysp), length(maxns))
grads     = zeros(length(xsp), length(ysp), length(maxns))
println("Looping through maxns...")
for (im,maxn) in enumerate(maxns)
    println("maxn = $maxn")
    tic = Base.time()
    interp = svd_IceSheetDEM.do_kriging(grid_output, geotable, varg; maxn)
    toc = Base.time() - tic
    m_interp = zeros(length(x),length(y))
    m_interp[ir_sim[ir_keep]] .= mean.(interp.Z)
    m_interps[:,:,im] = m_interp[xsp, ysp]
    grad = mgradient(m_interp)
    grads[:,:,im] = grad[xsp, ysp]
    ∑grad[im] = sum(grad)
    wallt[im] = toc
    tt = round(toc ./ 60, digits=1)
    println("Wall time: $tt minutes")
end
dict_to_save = (; maxns, wallt, m_interps, grads)
jldsave(joinpath(main_output_dir, "kriging_findmaxn.jld2"); dict_to_save...)
