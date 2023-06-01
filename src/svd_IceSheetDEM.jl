__precompile__(false)
module svd_IceSheetDEM

using JLD2, NetCDF, Glob, ProgressMeter, ArgParse

export read_model_data, parse_commandline, get_indices
export shortread, gdalwarp, save_netcdf, get_options

function parse_commandline(args)
    s = ArgParseSettings()
    @add_arg_table s begin
        "--λ", "--lambda"
            help     = "regularization parameter for the least squares problem"
            arg_type = Float64
            default  = 1e5
        "--r"
            help     = "truncation of SVD"
            arg_type = Int
            default  = 10
        "--res"
            help     = "resolution in m, currently available at '1200' or '1800'"
            arg_type = String
            default  = "1200"
        "--train_folder"
            help     = "folder where the training data 'usurf_*' is stored"
            arg_type = String
            default  = "data/"
        "--obs"
            help     = "file of observations that the SVD is fitted to"
            arg_type = String
            default  = "data/aerodem_g1200m_geoid_corrected_1978_1987_mean.nc"
        "--obs_band_name"
            help     = "name of the surface elevation band in the netcdf file"
            arg_type = String
            default  = "surface_altitude"
        "--save"
            help     = "save the output in an .nc file (option takes no argument)"
            action   = :store_true
    end
    return parse_args(args,s)
end

include("read_in.jl")
include("save_indices.jl")
include("gdal_helpers.jl")

end # module svd_IceSheetDEM
