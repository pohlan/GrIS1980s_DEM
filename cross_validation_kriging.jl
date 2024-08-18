using svd_IceSheetDEM, NetCDF, NCDatasets, Meshes, Statistics, StatsBase, GeoStats, CSV, DataFrames, Printf, JLD2, Distributions
import Plots

# set target directories
main_output_dir  = joinpath("output","geostats_interpolation")
fig_dir          = joinpath(main_output_dir, "figures")
mkpath(fig_dir)
atm_dh_dest_file = joinpath(dirname(obs_ATM_file), "grimp_minus_atm.csv")

outline_shp_file      = "data/gris-imbie-1980/gris-outline-imbie-1980_updated.shp"
# template_file = "data/training_data_it2_600/usurf_ex_gris_g600m_v2023_GIMP_id_0_1980-1-1_2020-1-1.nc"
template_file = "data/grimp/surface/grimp_geoid_corrected_g2400.nc"
gr = 2400

# get data
ref_file                    = "data/grimp/surface/grimp_geoid_corrected_g2400.nc"
reference_file_g150         = "data/grimp/surface/grimp_geoid_corrected_g150.nc"
bedm_file                   = create_bedmachine_grid(gr, template_file)
bedmachine_path             = splitdir(bedm_file)[1]
aerodem_g150, obs_aero_file = create_aerodem(;gr, outline_shp_file, bedmachine_path, reference_file_g150, kw="")
obs_ATM_file                = get_atm_file()
dhdt_file, _                = create_dhdt_grid(;gr, startyr=1994, endyr=2010)
mask_file                   = create_imbie_mask(;gr, outline_shp_file, sample_path=aerodem_g150)

# number of bins
nbins1 = 7
nbins2 = 9
blockspacing = gr*2

# preprocessing, standardization and variogram
df_all, varg, destand, I_no_ocean, idx_aero = svd_IceSheetDEM.prepare_interpolation(ref_file, bedm_file, obs_aero_file, obs_ATM_file, dhdt_file, mask_file, fig_dir, atm_dh_dest_file;
                                                                                    nbins1, nbins2, blockspacing)

geotable = svd_IceSheetDEM.make_geotable(df_all.dh_detrend, df_all.x, df_all.y)

xpl = [first(geom.coords) for geom in geotable.geometry]
ypl = [last(geom.coords) for geom in geotable.geometry]
Plots.scatter(xpl./1e3, ypl./1e3, marker_z=geotable.Z, label="", markersize=2.0, markerstrokewidth=0, cmap=:RdBu, clims=(-5,5), aspect_ratio=1, xlims=(-7e2,8e2), xlabel="Easting [km]", ylabel="Northing [km]", colorbar_title="[m]", title="dh standardized", grid=false, wsize=(1700,1800))

# create sets of training and test data
ℓ    = 2e4
flds = folds(geotable, BlockFolding(ℓ))
# flds = folds(geotable, BallFolding(MetricBall(ℓ)))
# flds = folds(geotable, OneFolding())

maxns = [300] #,800,1000]  #,2000,3000,4000,5000]

methods_name = ["median", "mean", "nmad", "std", "L2norm"]
methods_fct  = [median, mean, mad, std, norm]
dict = Dict{String,Array}(n => zeros(length(maxns)) for n in methods_name)
dict["maxn"] = maxns
m_difs = [Float32[] for i in eachindex(maxns)]
m_dists = [Float32[] for i in eachindex(maxns)]
for (im,maxn) in enumerate(maxns)
    println("maxn = $maxn")
    function evaluate_fun(i_train,i_test)
        interp, _ = svd_IceSheetDEM.do_kriging(view(domain(geotable),i_test), view(geotable,i_train), varg; maxn)
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
Plots.savefig(joinpath(fig_dir,"validation_histogram"))

# analyze error vs. distance to nearest neighbor, using binning
p = Plots.plot(xlabel="distance to nearest neighbor (km)", ylabel="mean error", legend_title="max # neighbors:", size=(1300,800), leftmargin=10Plots.mm, topmargin=10Plots.mm, bottommargin=10Plots.mm, legend=:top) #, palette = :Blues_7)
for (dists,difs,maxn) in zip(m_dists, m_difs, maxns) # eachrow(matrix_difs)
    dh_binned, bin_centers = svd_IceSheetDEM.bin_equal_sample_size(dists, difs, 4000)
    Plots.plot!(bin_centers[1:end-1]./1e3, median.(dh_binned[1:end-1]), marker=:circle, markersize=6, markerstrokewidth=0, label=maxn)
end
Plots.plot(p)
Plots.savefig(joinpath(fig_dir,"error_vs_distance.png"))
