using svd_IceSheetDEM, NetCDF, NCDatasets, Meshes, Statistics, StatsBase, GeoStats, CSV, DataFrames, Printf, JLD2, Distributions, ImageMorphology, StatsPlots, LaTeXStrings, Interpolations, UnPack
import Plots

# set target directories
main_output_dir = joinpath("output","validation")
fig_dir         = joinpath(main_output_dir, "figures")

outline_shp_file      = "data/gris-imbie-1980/gris-outline-imbie-1980_updated.shp"
grd = 600

################################################
# Preprocessing, standardization and variogram #
################################################

csv_preprocessing, jld2_preprocessing = svd_IceSheetDEM.prepare_obs(grd, outline_shp_file; blockspacing=grd/3, nbins1=7, nbins2=12)

# get I_no_ocean, (de-)standardization functions and variogram from pre-processing
df_all = CSV.read(csv_preprocessing, DataFrame)
dict   = load(jld2_preprocessing)
@unpack I_no_ocean, idx_aero, params = dict
@unpack destandardize = svd_IceSheetDEM.get_stddization_fcts(jld2_preprocessing)
varg = svd_IceSheetDEM.custom_var(params)

# make geotable
geotable = svd_IceSheetDEM.make_geotable(df_all.dh_detrend, df_all.x, df_all.y)

# xpl = [first(geom.coords) for geom in geotable.geometry]
# ypl = [last(geom.coords) for geom in geotable.geometry]
# Plots.scatter(xpl./1e3, ypl./1e3, marker_z=geotable.Z, label="", markersize=2.0, markerstrokewidth=0, cmap=:RdBu, clims=(-5,5), aspect_ratio=1, xlims=(-7e2,8e2), xlabel="Easting [km]", ylabel="Northing [km]", colorbar_title="[m]", title="dh standardized", grid=false, wsize=(1700,1800))

#######################################################
# Standard cross-validation, leave block/ball/one out #
#######################################################

# create sets of training and test data
ℓ    = 2e5
flds = folds(geotable, BlockFolding(ℓ))
# flds = folds(geotable, BallFolding(MetricBall(ℓ)))
# flds = folds(geotable, OneFolding())

maxn = 1600

function evaluate_fun(i_train,i_test)
    interp = svd_IceSheetDEM.do_kriging(view(domain(geotable),i_test), view(geotable,i_train), varg; maxn)
    return interp.Z
end
difs, xc, yc = svd_IceSheetDEM.step_through_folds(flds, evaluate_fun, geotable, save_coords=true, save_distances=false)

# get indices
idxs = [Int[] for i in flds]
for (i,fs) in enumerate(flds)
    idxs[i] = fs[2]
end
idx = vcat(idxs...)

logℓ      = round(log10(ℓ),digits=1)
dest = joinpath(main_output_dir,"cv_1e$(logℓ)_gr$(grd)_kriging.jld2")
jldsave(dest; difs, xc, yc, idx, binfield1=df_all.bfield_1[idx], h_ref=df_all.h_ref[idx], grd, method="kriging")



#######################################################

# plot map of difs
p = Plots.scatter(xc,yc,marker_z=difs, cmap=:bwr, clims=(-4,4), markersize=0.7, markerstrokewidth=0, label="", aspect_ratio=1, size=(500,700))
# p = Plots.plot()
for (j,fs) in enumerate(flds)
    stest = view(domain(geotable), fs[2])
    ch    = GeoStats.convexhull(stest).rings[1].vertices.data
    xx    = first.(a.coords for a in ch)
    yy    = last.(a.coords for a in ch)
    Plots.plot!(p, xx,yy, lw=1, color="black", label="")
end
Plots.plot(p, size=(500,700))
Plots.savefig(joinpath(fig_dir, "map_validation_maxn$maxn.png"))

#######################################################
# Determine good maxns through morphological gradient #
#######################################################
bedmachine_original, bedm_file = svd_IceSheetDEM.create_bedmachine_grid(grd)
reference_file_g150, ref_file  = svd_IceSheetDEM.create_grimpv2(grd, bedmachine_original)

# derive indices for cells to interpolate, this time interpolate the points we don't know, technically no validation!
ir_sim      = setdiff(I_no_ocean, idx_aero)  # indices that are in I_no_ocean but not in idx_aero
xsp = 800:1150
ysp = 500:850
x           = NCDataset(ref_file)["x"][:]
y           = NCDataset(ref_file)["y"][:]
ix = get_ix.(ir_sim, length(x))
iy = get_iy.(ir_sim, length(x))
ir_keep = findall(xsp[1].<ix.<xsp[end] .&& ysp[1].<iy.<ysp[end])
grid_output = PointSet([Point(xi,yi) for (xi,yi) in zip(x[ix[ir_keep]], y[iy[ir_keep]])])

ir_keep_df = findall(x[xsp[1]].<df_all.x.<x[xsp[end]] .&& y[ysp[1]].<df_all.y.<y[ysp[end]])
geotable = svd_IceSheetDEM.make_geotable(df_all.dh_detrend[ir_keep_df], df_all.x[ir_keep_df], df_all.y[ir_keep_df])


# loop through maxns and calculate the morphological gradient each time; if high, good indication that maxn is too low
maxns = [1500, 2000, 2300]
∑grad = zeros(length(maxns))
wallt = zeros(length(maxns))
for (im,maxn) in enumerate(maxns)
    println("maxn = $maxn")
    tic = Base.time()
    interp = svd_IceSheetDEM.do_kriging(grid_output, geotable, varg; maxn)
    toc = Base.time() - tic
    m_interp = zeros(length(x),length(y))
    m_interp[ir_sim[ir_keep]] .= mean.(interp.Z)
    grad = mgradient(m_interp)
    ∑grad[im] = sum(grad)
    wallt[im] = toc
    p1 = Plots.heatmap(m_interp[xsp,ysp]',cmap=:bwr, clims=(-4,4), aspect_ratio=1)
    p2 = Plots.heatmap(grad[xsp,ysp]', aspect_ratio=1)
    Plots.plot(p1,p2)
    Plots.savefig(joinpath(fig_dir, "kriging_maxn$(maxn)_g$(grd).png"))
    tt = round(toc ./ 60, digits=1)
    sgrd = ∑grad[im]
    print("Wall time: $tt minutes, "); println("sumgrad = $sgrd")
end
Plots.plot(maxns, ∑grad, marker=:circle, label="", xlabel="max  # neighbors kriging", title="Beucher morphological gradient", ylabel="sum(gradient)")
Plots.savefig(joinpath(fig_dir, "maxn_vs_gradient.png"))
