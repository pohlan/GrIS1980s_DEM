using svd_IceSheetDEM, NetCDF, NCDatasets, Meshes, Statistics, StatsBase, GeoStats, CSV, DataFrames, Printf, JLD2, Distributions, ImageMorphology, StatsPlots, LaTeXStrings, Interpolations, UnPack
import Plots

# set target directories
main_output_dir = joinpath("output","validation")

outline_shp_file      = "data/gris-imbie-1980/gris-outline-imbie-1980_updated.shp"
gr = 1200

################################################
# Preprocessing, standardization and variogram #
################################################

csv_preprocessing, jld2_preprocessing = svd_IceSheetDEM.prepare_obs(gr, outline_shp_file; blockspacing=gr/3, nbins1=7, nbins2=12)

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

maxn = 100

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
dest = joinpath(main_output_dir,"dict_cv_block_1e$(logℓ)_gr$(gr)_kriging.jld2")
jldsave(dest; difs, xc, yc, idx, h_ref=df_all.h_ref[idx], gr, method="kriging")



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

# derive indices for cells to interpolate, this time interpolate the points we don't know, technically no validation!
ir_sim      = setdiff(I_no_ocean, idx_aero)  # indices that are in I_no_ocean but not in idx_aero
x           = NCDataset(obs_aero_file)["x"][:]
y           = NCDataset(obs_aero_file)["y"][:]
grid_output = PointSet([Point(xi,yi) for (xi,yi) in zip(x[get_ix.(ir_sim, length(x))], y[get_iy.(ir_sim, length(x))])])
geotable    = svd_IceSheetDEM.make_geotable(df_all.dh_detrend, df_all.x, df_all.y)

# loop through maxns and calculate the morphological gradient each time; if high, good indication that maxn is too low
maxns = [10,50,100,200,300,500]
∑grad = zeros(length(maxns))
for (im,maxn) in enumerate(maxns)
    tic = Base.time()
    interp = svd_IceSheetDEM.do_kriging(grid_output, geotable, varg; maxn)
    toc = Base.time() - tic
    m_interp = zeros(length(x),length(y))
    m_interp[ir_sim] .= interp.Z
    grad = mgradient(m_interp)
    ∑grad[im] = sum(grad)
    p1 = Plots.heatmap(m_interp[200:600,1:600]',cmap=:bwr, clims=(-4,4))
    p2 = Plots.heatmap(grad[200:600,1:600]')
    Plots.plot(p1,p2)
    Plots.savefig(joinpath(fig_dir_krig, "kriging_maxn$(maxn)_g(gr).png"))
end
Plots.plot(maxns, ∑grad, marker=:circle, label="", xlabel="max  # neighbors kriging", title="Beucher morphological gradient", ylabel="sum(gradient)")
Plots.savefig(joinpath(fig_dir_krig, "maxn_vs_gradient.png"))
