__precompile__(false)
module svd_IceSheetDEM

using ArgParse

export parse_commandline
export shortread, gdalwarp
export create_aerodem, create_bedmachine_grid, create_imbie_mask
export solve_lsqfit, create_reconstructed_bedmachine

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
        "--train_folder"
            help     = "folder where the training data 'usurf_*' is stored"
            arg_type = String
            default  = "data/training_data_it0_1200/"
        "--imbie_path"
            help     = "folder where imbie shp file is stored"
            arg_type = String
            default  = "data/gris-imbie-1980/"
    end
    return parse_args(args,s)
end

include("read_in.jl")
include("gdal_helpers.jl")
include("rsvd.jl")
include("reconstruction_routines.jl")

end # module svd_IceSheetDEM
