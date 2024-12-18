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

csv_preprocessing, jld2_preprocessing = svd_IceSheetDEM.prepare_obs(grd, outline_shp_file)

# get I_no_ocean, (de-)standardization functions and variogram from pre-processing
df_all = CSV.read(csv_preprocessing, DataFrame)
dict   = load(jld2_preprocessing)
@unpack I_no_ocean, idx_aero, gamma = dict
@unpack destandardize = svd_IceSheetDEM.get_stddization_fcts(jld2_preprocessing)
# varg = svd_IceSheetDEM.custom_var(params)
varg = svd_IceSheetDEM.get_var(gamma)

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

maxns = [300]

for maxn in maxns
    println(maxn)
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
    dest = joinpath(main_output_dir,"cv_1e$(logℓ)_gr$(grd)_kriging_maxn$(maxn).jld2")
    jldsave(dest; maxn, difs, xc, yc, idx, binfield1=df_all.bfield_1[idx], h_ref=df_all.h_ref[idx], grd, method="kriging")
end



#######################################################

# plot map of difs
p = Plots.scatter(xc,yc,marker_z=difs, cmap=:bwr, clims=(-4,4), markersize=0.7, markerstrokewidth=0, label="", aspect_ratio=1, size=(500,700))
# p = Plots.plot()
for (j,fs) in enumerate(flds)
    stest = view(domain(geotable), fs[2])
    ch    = GeoStats.convexhull(stest).rings[1].vertices.data
    xx    = [ustrip.(a.coords.x) for a in ch]
    yy    = [ustrip.(a.coords.y) for a in ch]
    Plots.plot!(p, xx,yy, lw=1, color="black", label="")
end
Plots.plot(p, size=(500,700))
Plots.savefig(joinpath(fig_dir, "map_validation_maxn$maxn.png"))

#######################################################
# Determine good maxns through morphological gradient #
#######################################################
bedmachine_original, bedm_file = svd_IceSheetDEM.create_bedmachine_grid(grd)
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


# loop through maxns and calculate the morphological gradient each time; if high, good indication that maxn is too low
maxns = [10, 100, 1500, 3000]∑grad = zeros(length(maxns))
wallt = zeros(length(maxns))
m_interps = zeros(length(xsp), length(ysp), length(maxns))
grads     = zeros(length(xsp), length(ysp), length(maxns))
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
    p1 = Plots.heatmap(m_interp[xsp,ysp]',cmap=:bwr, clims=(-4,4), aspect_ratio=1)
    p2 = Plots.heatmap(grad[xsp,ysp]', aspect_ratio=1)
    Plots.plot(p1,p2)
    Plots.savefig(joinpath(fig_dir, "kriging_maxn$(maxn)_g$(grd).png"))
    tt = round(toc ./ 60, digits=1)
    sgrd = ∑grad[im]
    print("Wall time: $tt minutes, "); println("sumgrad = $sgrd")
end
dict_to_save = (; maxns, wallt, m_interps, grads)
jldsave("output/validation/kriging_findmaxn.jld2"; dict_to_save...)

Plots.plot(maxns, ∑grad, marker=:circle, label="", xlabel="max  # neighbors kriging", title="Beucher morphological gradient", ylabel="sum(gradient)")
Plots.savefig(joinpath(fig_dir, "maxn_vs_gradient.png"))


outline_shp_file = "data/gris-imbie-1980/gris-outline-imbie-1980_updated_crs.shp"
shp              = Shapefile.shapes(Shapefile.Table(outline_shp_file))

Plots.scalefontsizes()
Plots.scalefontsizes(1.6)
@unpack maxns, m_interps = load("output/validation/kriging_findmaxn.jld2")
ps = Plots.Plot{Plots.GRBackend}[]
clims=(-4,4)
cmap = :RdBu
for m in axes(m_interps, 3)
    # cbar = m == 3 ? true : false
    yticks = m == 1 ? true : false
    pi = heatmap(x[xsp], y[ysp], m_interps[:,:,m]', title="\n"*L"\mathrm{n_{obs}} = "*"$(maxns[m])", aspect_ratio=1, wsize=(1500,400), grid=false, cbar=false, tick_direction=:out, titlefontsize=18; clims, cmap)
    if m !== 1
        pi = plot(pi, ytickfontsize=1, ytickfontcolor=:white)
    end
    plot!(shp, xlims=extrema(x[xsp]), ylims=extrema(y[ysp]), fill=nothing, lw=0.5)
    push!(ps, pi)
end
p_panels = plot(ps..., size=(3000, 500), layout=(1,4), leftmargin=10Plots.mm, rightmargin=10Plots.mm, topmargin=-10Plots.mm, bottommargin=-10Plots.mm)

xx = range(clims...,1000)
zz = zero(xx)' .+ xx
p_c = heatmap(xx, xx, zz, ratio=15, xticks=false, legend=false, fc=cgrad(cmap), lims=(-4,4), framestyle=:box, left_margin=-200Plots.mm, top_margin=30Plots.mm, bottom_margin=30Plots.mm, ymirror=true) #, size=(10,100))
annotate!(40, 0.0, text("m", 18, "Computer Modern", :left))
plot(p_c)

# plot again everything together
plot(p_panels, p_c, right_margin=-10Plots.mm) #; bottom_margin=-40Plots.mm, size=(2100,600), top_margin=10Plots.mm)
savefig("output/validation/figures/kriging_interp_maps_maxn.png")
