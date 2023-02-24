__precompile__(false)
module svd_IceSheetDEM

using JLD2, NetCDF, Glob, ProgressMeter

export read_model_data

include("read_in.jl")

if !all(isfile.(["data/indices_aerodem_g1200m.jld2", "data/indices_aerodem_g1800m.jld2"]))
    include("save_indices.jl")
end

end # module svd_IceSheetDEM
