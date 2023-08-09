__precompile__(false)
module svd_IceSheetDEM

using ArgParse

export parse_commandline
export archgdal_read, gdalwarp
export create_aerodem, create_bedmachine_grid, create_imbie_mask
export solve_lsqfit, create_reconstructed_bedmachine
export rsvd

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
            required = true
        "--imbie_shp_file"
            help     = "shape file outlining the ice"
            arg_type = String
    end
    return parse_args(args,s)
end

include("read_in.jl")
include("gdal_helpers.jl")
include("rsvd.jl")
include("reconstruction_routines.jl")

end # module svd_IceSheetDEM
