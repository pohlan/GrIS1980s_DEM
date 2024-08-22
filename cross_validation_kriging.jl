using svd_IceSheetDEM, NetCDF, NCDatasets, Meshes, Statistics, StatsBase, GeoStats, CSV, DataFrames, Printf, JLD2, Distributions, ImageMorphology
import Plots

# set target directories
main_output_dir  = joinpath("output","geostats_interpolation")
fig_dir          = joinpath(main_output_dir, "figures")
fig_dir_krig     = joinpath(main_output_dir, "kriging", "figures")
mkpath(fig_dir)

outline_shp_file      = "data/gris-imbie-1980/gris-outline-imbie-1980_updated.shp"
gr = 1200

# get filenames
bedmachine_original, bedm_file = create_bedmachine_grid(gr)
reference_file_g150, ref_file  = create_grimpv2(gr, bedmachine_original)
aerodem_g150, obs_aero_file    = create_aerodem(gr, outline_shp_file, bedmachine_original, reference_file_g150)
obs_ATM_file                   = get_atm_file()
dhdt_file, _                   = create_dhdt_grid(;gr, startyr=1994, endyr=2010)
mask_file                      = create_imbie_mask(;gr, outline_shp_file, sample_path=aerodem_g150)
atm_dh_dest_file               = joinpath(dirname(obs_ATM_file), "grimp_minus_atm.csv")

################################################
# Preprocessing, standardization and variogram #
################################################

# number of bins
nbins1 = 7
nbins2 = 12
blockspacing = gr

df_all, varg, destand, I_no_ocean, idx_aero = svd_IceSheetDEM.prepare_interpolation(ref_file, bedm_file, obs_aero_file, obs_ATM_file, dhdt_file, mask_file, fig_dir, atm_dh_dest_file;
                                                                                    nbins1, nbins2, blockspacing)

geotable = svd_IceSheetDEM.make_geotable(df_all.dh_detrend, df_all.x, df_all.y)
xpl = [first(geom.coords) for geom in geotable.geometry]
ypl = [last(geom.coords) for geom in geotable.geometry]
Plots.scatter(xpl./1e3, ypl./1e3, marker_z=geotable.Z, label="", markersize=2.0, markerstrokewidth=0, cmap=:RdBu, clims=(-5,5), aspect_ratio=1, xlims=(-7e2,8e2), xlabel="Easting [km]", ylabel="Northing [km]", colorbar_title="[m]", title="dh standardized", grid=false, wsize=(1700,1800))

#######################################################
# Standard cross-validation, leave block/ball/one out #
#######################################################

# create sets of training and test data
ℓ    = 2e4
flds = folds(geotable, BlockFolding(ℓ))
# flds = folds(geotable, BallFolding(MetricBall(ℓ)))
# flds = folds(geotable, OneFolding())

maxns = [10,50,100,500,800,1000]

methods_name = ["median", "mean", "nmad", "std", "L2norm"]
methods_fct  = [median, mean, mad, std, norm]
dict = Dict{String,Array}(n => zeros(length(maxns)) for n in methods_name)
dict["maxn"] = maxns
m_difs = [Float32[] for i in eachindex(maxns)]
m_dists = [Float32[] for i in eachindex(maxns)]
for (im,maxn) in enumerate(maxns)
    println("maxn = $maxn")
    function evaluate_fun(i_train,i_test)
        interp = svd_IceSheetDEM.do_kriging(view(domain(geotable),i_test), view(geotable,i_train), varg; maxn)
        return interp.Z
    end
    tic = Base.time()
    difs, dists = svd_IceSheetDEM.step_through_folds(flds, evaluate_fun, geotable, save_distances=true)
    toc = Base.time() - tic
    @printf("Kriging took %1.1f minutes. \n", toc / 60)
    for (mn, mf) in zip(methods_name, methods_fct)
        dict[mn][im] = mf(difs)
    end
    # save dists and difs
    m_difs[im]  = difs
    m_dists[im] = dists
end
logℓ = round(log(10,ℓ),digits=1)
dest      = joinpath(main_output_dir,"dict_cv_block_$(logℓ).jld2")
save(dest, dict)

# plot histogram
ph = Plots.plot()
for difs in m_difs
    Plots.histogram!(ph, difs, normalize=:pdf, linecolor=nothing)
end
Plots.savefig(joinpath(fig_dir_krig,"validation_histogram"))

# analyze error vs. distance to nearest neighbor, using binning
p = Plots.plot(xlabel="distance to nearest neighbor (km)", ylabel="mean error", legend_title="max # neighbors:", size=(1300,800), leftmargin=10Plots.mm, topmargin=10Plots.mm, bottommargin=10Plots.mm, legend=:top) #, palette = :Blues_7)
for (dists,difs,maxn) in zip(m_dists, m_difs, maxns) # eachrow(matrix_difs)
    dh_binned, bin_centers = svd_IceSheetDEM.bin_equal_sample_size(dists, difs, 4000)
    Plots.plot!(bin_centers[1:end-1]./1e3, median.(dh_binned[1:end-1]), marker=:circle, markersize=6, markerstrokewidth=0, label=maxn)
end
Plots.plot(p)
Plots.savefig(joinpath(fig_dir_krig,"error_vs_distance.png"))

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
