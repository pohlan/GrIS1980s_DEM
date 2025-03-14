using svd_IceSheetDEM, NCDatasets, JLD2, UnPack, CSV, DataFrames, Glob, Dates, GeoDataFrames, GeoFormatTypes
using Plots, StatsPlots, LaTeXStrings, GeoStats, Shapefile, StatsBase, Meshes, Distributions

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

# load data
csv_preprocessing, jld2_preprocessing = svd_IceSheetDEM.prepare_obs(grd, "")
@unpack standardize, destandardize = svd_IceSheetDEM.get_stddization_fcts(jld2_preprocessing)
@unpack href_file, bin_centers_1, bin_centers_2, nmads, meds, gamma, I_no_ocean = load(jld2_preprocessing)
df_all = CSV.read(csv_preprocessing, DataFrame)

h_ref = NCDataset(href_file)["surface"][:]
x     = NCDataset(href_file)["x"][:]
y     = NCDataset(href_file)["y"][:]

# to plot outline, polygon needs to be reprojected
shp             = Shapefile.shapes(Shapefile.Table(outline_shp_file))
coords_inverted = [(pt.y, pt.x) for pt in shp[1].points]   # reproject of GeoDataFrames expects the coordinates in reverse order
df_inverted     = DataFrame(geometry=AG.createpolygon(coords_inverted))
reproject!(df_inverted, EPSG(4326), EPSG(3413))
outl            = df_inverted.geometry

include("plot_preprocessing.jl")
include("plot_eigenmodes.jl")
include("plot_validation_results.jl")
include("plot_glacier_flowlines.jl")
include("plot_initialization.jl")
