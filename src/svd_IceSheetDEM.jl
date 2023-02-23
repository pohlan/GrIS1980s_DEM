__precompile__(false)
module svd_IceSheetDEM

using JLD2, NetCDF, Glob

export read_model_data

include("read_in.jl")

if !isfile("data/aerodem_data_indices.jld2")
    include("save_indices.jl")
end

end # module svd_IceSheetDEM
