using svd_IceSheetDEM, NCDatasets, JLD2, UnPack, CSV, DataFrames, Glob, Dates, GeoDataFrames
using Plots, StatsPlots, LaTeXStrings, GeoStats, Shapefile, StatsBase, Statistics, Meshes, Distributions
import ArchGDAL as AG
import GeoFormatTypes as GFT

# set target directories for paper figures
fig_dir_main = joinpath("output","main_figures")
mkpath(fig_dir_main)

# for running the script interactively
# ARGS = [
#         "--shp_file", "data/gris-imbie-1980/gris-outline-imbie-1980_updated.shp",
#         "--grid_size", "600.0"]

parsed_args         = parse_commandline(ARGS)
outline_shp_file    = parsed_args["shp_file"]
grd                 = parsed_args["grid_size"]
init_run_krig       = parsed_args["initialization_run_kriging"]
init_run_SVD        = parsed_args["initialization_run_SVD"]

# load data
csv_preprocessing, jld2_preprocessing = svd_IceSheetDEM.prepare_obs(grd, "")
@unpack standardize, destandardize = svd_IceSheetDEM.get_stddization_fcts(jld2_preprocessing)
@unpack href_file, bin_centers_1, bin_centers_2, nmads, meds, gamma, I_no_ocean = load(jld2_preprocessing)
df_all = CSV.read(csv_preprocessing, DataFrame)

h_ref_m = nomissing(NCDataset(href_file)["surface"][:,:], 0.0)
x       = NCDataset(href_file)["x"][:]
y       = NCDataset(href_file)["y"][:]

# to plot outline, polygon needs to be reprojected
shp    = Shapefile.shapes(Shapefile.Table(outline_shp_file))
coords = [(pt.x, pt.y) for pt in shp[1].points]
df     = DataFrame(geometry=AG.createpolygon(coords))
reproject!(df, GFT.EPSG(4326), GFT.EPSG(3413), always_xy=true)
outl   = df.geometry

# chosen parameters for kriging and SVD
const maxn0 = 1500
const r0    = 500
const λ0    = 1e7
const logλ = Int(round(log(10, λ0)))

include("plot_preprocessing.jl")
include("plot_eigenmodes.jl")
include("plot_validation_results.jl")
include("plot_glacier_flowlines.jl")
include("plot_initialization.jl")
