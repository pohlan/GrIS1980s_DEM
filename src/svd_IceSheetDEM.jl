__precompile__(false)
module svd_IceSheetDEM

using ArgParse

import ArchGDAL as AG
import GeoFormatTypes as GFT

using DelimitedFiles, NCDatasets, NetCDF, Glob, DataFrames, CSV, Dates, GeoFormatTypes, ZipFile, JLD2, UnPack
using Downloads, Cascadia, Gumbo, HTTP, PyCall
using Printf, ProgressMeter
using Statistics, GeoStats, StatsBase, Distributions, Interpolations, LsqFit, ImageFiltering, ParallelRandomFields.grf2D_CUDA
using Arpack, LinearAlgebra, ImageMorphology
using DataStructures: OrderedDict
import Plots, StatsPlots

export parse_commandline, get_ix, get_iy
export SVD_reconstruction, create_reconstructed_bedmachine, prepare_obs
export SVD_random_fields, geostats_interpolation

const no_data_value = -9999.0
const F = Float32              # Julia default is Float64 but that kills the process for the full training data set if r is too large

function parse_commandline(args)
    s = ArgParseSettings()
    @add_arg_table s begin
        "--λ", "--lambda"
            help     = "regularization parameter for the least squares problem"
            arg_type = Float64
            default  = 1e5
        "--r"
            help     = "truncation of SVD, default is a full SVD"
            arg_type = Int
            default  = 10^30  # numbers close to the maximum or larger will give a full SVD
        "--training_data"
            help     = "training files, e.g. train_folder/usurf*.nc"
            nargs    = '*'
            arg_type = String
            required = true
        "--shp_file"
            help     = "shape file outlining the ice"
            arg_type = String
            required = true
        "--do_figures"
            help     = "whether or not to plot the difference of the reconstructed elevations to aerodem data (plotting and saving requires a bit of extra memory)"
            arg_type = Bool
            default  = false
        "--use_arpack"
            help     = "If this is set to true, then the Arpack svd is used instead of the standard LinearAlgebra algorithm. Arpack is iterative and matrix free and thus useful when memory becomes limiting, but it can be slower."
            arg_type = Bool
            default  = false
        "--maxn"
            help     = "Maximum number of neighbors used for kriging."
            arg_type = Int
            default  = 1600
    end
    return parse_args(args,s)
end

get_ix(i,nx) = i % nx == 0 ? nx : i % nx
get_iy(i,nx) = cld(i,nx)
get_global_i(ix, iy, nx) = nx * (iy-1) + ix

include("prepare_data.jl")
include("gdal_helpers.jl")
include("reconstruction_routines.jl")
include("statistics_helpers.jl")
include("plotting_tools.jl")

end # module svd_IceSheetDEM
