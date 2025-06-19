__precompile__(false)
module GrIS1980s_DEM

using ArgParse

import ArchGDAL as AG
import GeoFormatTypes as GFT

using DelimitedFiles, NCDatasets, NetCDF, Glob, DataFrames, CSV, Dates, GeoFormatTypes, ZipFile, JLD2, UnPack
using Downloads, Cascadia, Gumbo, HTTP, PythonCall
using Printf, ProgressMeter
using Statistics, GeoStats, StatsBase, Distributions, Interpolations, ImageFiltering, LocalFilters
using Arpack, LinearAlgebra, ImageMorphology, Optim
using DataStructures: OrderedDict
using Plots, StatsPlots, LaTeXStrings

export parse_commandline, get_ix, get_iy
export get_std_uncrt_file, get_cv_file_SVD, get_cv_file_kriging, get_rec_file_SVD, get_rec_file_kriging, get_rec_file_SVD_combined, kriging_findmaxn_file
export uncertainty_from_cv, create_reconstructed_bedmachine, prepare_obs, download_velocity
export SVD_reconstruction, geostats_interpolation, combined_SVD_AeroDEM

const no_data_value = -9999.0
const F = Float32              # Julia default is Float64 but that kills the process for the full training data set if r is too large

function parse_commandline(args)
    s = ArgParseSettings()
    @add_arg_table s begin
        "--λ", "--lambda"
            help     = "regularization parameter for the least squares problem"
            arg_type = Float64
            default  = 1e7
        "--r"
            help     = "truncation of SVD, default is a full SVD"
            arg_type = Int
            default  = 10^30  # numbers close to the maximum or larger will give a full SVD
        "--training_data"
            help     = "model-generated realizations of ice sheet elevation as netcdf files, e.g. train_folder/usurf*.nc"
            nargs    = '*'
            arg_type = String
        "--initialization_run_kriging"
            help     = "simulation initialized with kriging DEM"
            arg_type = String
        "--initialization_run_SVD"
            help     = "simulation initialized with SVD DEM"
            arg_type = String
        "--shp_file"
            help     = "shape file outlining the ice sheet"
            arg_type = String
        "--use_arpack"
            help     = "if set to true, the Arpack svd is used instead of the standard LinearAlgebra algorithm; Arpack is iterative and matrix free and thus useful when memory becomes limiting, but it can be slower"
            arg_type = Bool
            default  = false
        "--maxn"
            help     = "maximum number of neighbors used for kriging"
            arg_type = Int
            default  = 1500
        "--grid_size"
            help     = "cell size of grid, same in x and y direction; not needed for svd reconstruction where training_data is provided"
            arg_type = Int
            default  = 600.0
        "--svd_mask"
            help     = "for combined SVD/AeroDEM: shape file specifying where SVD rec. should be used"
            arg_type = String
    end
    return parse_args(args,s)
end

get_ix(i,nx) = i % nx == 0 ? nx : i % nx
get_iy(i,nx) = cld(i,nx)
get_global_i(ix, iy, nx) = nx * (iy-1) + ix

include("gdal_helpers.jl")
include("statistics_helpers.jl")
include("kriging_interpolation.jl")
include("SVD_reconstruction.jl")
include("plotting_tools.jl")

end # module GrIS1980s_DEM
