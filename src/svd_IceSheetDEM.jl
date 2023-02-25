__precompile__(false)
module svd_IceSheetDEM

using JLD2, NetCDF, Glob, ProgressMeter, ArgParse

export read_model_data, parse_commandline

function parse_commandline()
    s = ArgParseSettings()
    @add_arg_table s begin
        "--λ", "--lambda"
            help     = "regularization parameter for the least squares problem"
            arg_type = Float64
            default  = 0.
        "--r", "--trunc"
            help     = "number of modes that should be kept in the SVD truncation"
            arg_type = Int
            default  = nothing
        "--res"
            help     = "fesolution in m, currently available at '1200' or '1800'"
            arg_type = String
            default  = "1200"
        "--filepath"
            help     = "folder where the training data 'usurf_*' is stored"
            arg_type = String
            default  = "data/"
        "--save"
            help     = "save the output in an .nc file (option takes no argument)"
            action   = :store_true
        # doesn't currently work if one wants to pass an actual argument, for --F "Float32" it gives the error
        # invalid argument: Float32 (conversion to type DataType failed; you may need to overload ArgParse.parse_item;
        # "--F"
        #     help     = "Precision, Julia default is Float64 but that kills the process for the full training data set if r is too large"
        #     arg_type = DataType
        #     default  = Float32
    end
    return parse_args(s)
end

include("read_in.jl")

if !all(isfile.(["data/indices_aerodem_g1200m.jld2", "data/indices_aerodem_g1800m.jld2"]))
    include("save_indices.jl")
end

end # module svd_IceSheetDEM
