using GrIS1980s_DEM, NetCDF, NCDatasets, Meshes, Statistics, StatsBase, GeoStats, CSV, DataFrames, Printf, JLD2, Distributions, ImageMorphology, StatsPlots, LaTeXStrings, Interpolations, UnPack
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

################################################
# Preprocessing, standardization and variogram #
################################################

csv_preprocessing, jld2_preprocessing = GrIS1980s_DEM.prepare_obs(grd, outline_shp_file)

# get I_no_ocean, (de-)standardization functions and variogram from pre-processing
df_all = CSV.read(csv_preprocessing, DataFrame)
dict   = load(jld2_preprocessing)
@unpack I_no_ocean, idx_aero, gamma, href_file = dict
@unpack destandardize = GrIS1980s_DEM.get_stddization_fcts(jld2_preprocessing)
varg = GrIS1980s_DEM.get_var(gamma)

# make geotable
geotable = GrIS1980s_DEM.make_geotable(df_all.dh_detrend, df_all.x, df_all.y)

##############################################
# Standard cross-validation, leave block out #
##############################################

# create sets of training and test data
ℓ    = 2e5
flds = folds(geotable, BlockFolding(ℓ))

maxns = [10, 100, 500, 1500, 3000]
println("Kriging cross-validation...")
for maxn in maxns
    println("maxn = $(maxn)...")
    function evaluate_fun(i_train,i_test)
        interp = GrIS1980s_DEM.do_kriging(view(domain(geotable),i_test), view(geotable,i_train), varg; maxn)
        return interp.Z
    end
    difs, xc, yc = GrIS1980s_DEM.step_through_folds(flds, evaluate_fun, geotable, save_coords=true, save_distances=false)

    # get indices
    idxs = [Int[] for i in flds]
    for (i,fs) in enumerate(flds)
        idxs[i] = fs[2]
    end
    idx = vcat(idxs...)

    logℓ    = round(log10(ℓ),digits=1)
    dest    = get_cv_file_kriging(grd, maxn, logℓ)
    cv_dict = (; maxn, difs, xc, yc, idx, binfield1=df_all.bfield_1[idx], h_ref=df_all.h_ref[idx], grd, method="kriging")
    jldsave(dest; cv_dict...)

    # plot map of difs and blocks
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
end


###############################################
# Interpolate a subregion for different maxns #
###############################################
bedmachine_original, bedm_file    = GrIS1980s_DEM.create_bedmachine_grid(grd)
reference_file_g150, _, ref_file  = GrIS1980s_DEM.create_grimpv2(grd, dict["coreg_grid"], bedmachine_original)

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
geotable = GrIS1980s_DEM.make_geotable(df_all.dh_detrend[ir_keep_df], df_all.x[ir_keep_df], df_all.y[ir_keep_df])

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
    interp = GrIS1980s_DEM.do_kriging(grid_output, geotable, varg; maxn)
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
dict_to_save = (; xsp, ysp, maxns, wallt, m_interps, grads)
jldsave(kriging_findmaxn_file(); dict_to_save...)
