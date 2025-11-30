using GrIS1980s_DEM, Glob, NCDatasets, JLD2, GeoStats, Statistics, StatsBase, StatsPlots, LaTeXStrings, CSV, DataFrames, UnPack, LinearAlgebra
import Plots

# set target directories
main_output_dir = joinpath("output","validation")
mkpath(main_output_dir)

# for running the script interactively
# ARGS = [
#         "--shp_file", "data/gris-imbie-1980/gris-outline-imbie-1980_updated.shp",
#         "--training_data", readdir("data/training_data_it2_600", join=true)...]

parsed_args = parse_commandline(ARGS)
model_files = parsed_args["training_data"]
outline_shp_file = parsed_args["shp_file"]
use_arpack  = parsed_args["use_arpack"]

only_atm = false

# determine grid resolution
template_file = model_files[1]
x = NCDataset(template_file)["x"][:]
y = NCDataset(template_file)["y"][:]
const grd = Int(x[2] - x[1])

csv_preprocessing, jld2_preprocessing = GrIS1980s_DEM.prepare_obs(grd, outline_shp_file)

# get I_no_ocean, (de-)standardization functions and variogram from pre-processing
df_all = CSV.read(csv_preprocessing, DataFrame)
dict   = load(jld2_preprocessing)
@unpack I_no_ocean, idx_aero, href_file = dict

# load datasets, take full SVD (to be truncated with different rs later) for different numbers of training data files
# This step saves the output to files and is skipped if the files exist already; delete files to force
n_training_files = [10,30,50,70]
fnames_modes = ["" for i in eachindex(n_training_files)]
for (i,nfils) in enumerate(n_training_files)
    fname = joinpath(main_output_dir, "SVD_components_g$(grd)_nfiles$(nfils).jld2")
    fnames_modes[i] = fname
    if isfile(fname)
        continue
    end
	GrIS1980s_DEM.prepare_model(model_files[1:nfils], I_no_ocean, fname) # read in model data and take svd to derive "eigen ice sheets
end

# get I_obs and i_atm (will be the same for all nfiles)
@unpack data_mean = load(fnames_modes[1])
x_data, I_obs, i_atm = GrIS1980s_DEM.prepare_obs_SVD(grd, csv_preprocessing, I_no_ocean, data_mean, main_output_dir)

# create geotable (for partitioning)
x_Iobs   = x[get_ix.(I_no_ocean[I_obs],length(x))]
y_Iobs   = y[get_iy.(I_no_ocean[I_obs],length(x))]
coords_Iobs = [[ix,iy] for (ix,iy) in zip(x_Iobs,y_Iobs)]
geotable = georef(nothing, coords_Iobs)
geotable_atm = georef(nothing, coords_Iobs[i_atm])

# create sets of training and test data
ℓ    = 2e5 # radius for leaving out Block
logℓ = round(log(10,ℓ),digits=1)
geotable_partition = only_atm ? geotable_atm : geotable
blocks = partition(geotable_partition, BlockPartition(ℓ))
ids_test  = indices(blocks)
ids_train = []
for i in eachindex(ids_test)
    if only_atm ids_test[i] = i_atm[ids_test[i]] end
    id_train    = Vector(1:length(I_obs))
    deleteat!(id_train, ids_test[i])
    push!(ids_train, id_train)
end

# calculate distances to nearest observation
dists = get_closest_neighbors_from_cv(ids_train, ids_test, x_Iobs, y_Iobs)

# give λ and r values to loop through
λs        = [1e5, 1e6, 1e7, 1e8, 1e9, 1e10, 1e11]
rs        = [10, 50, 100, 300, 500, 700]
function do_validation_and_save(f)
    # load data
    @unpack U, Σ, data_mean, nfiles = load(f)
    UΣ = U*diagm(Σ)

    # load datasets, take full SVD (to be truncated with different rs later)
    x_data, _, _ = GrIS1980s_DEM.prepare_obs_SVD(grd, csv_preprocessing, I_no_ocean, data_mean, main_output_dir)

    function predict_vals(λ, r, i_train, i_test, x_data, I_obs, UΣ)
        _, x_rec = GrIS1980s_DEM.solve_optim(UΣ, I_obs[i_train], r, λ, x_data[i_train])
        return x_rec[I_obs[i_test]]
    end

    # loop through λ and r values
    m_difs  = [Float32[] for i in eachindex(λs), j in eachindex(rs)]
    for (iλ,λ) in enumerate(λs)
        for (ir,r) in enumerate(rs)
            if r > size(UΣ,2) continue end
            logλ = round(log(10, λ),digits=1)
            println("r = $r, logλ = $logλ")
            evaluate_fun(i_train,i_test) = predict_vals(λ, r, i_train, i_test, x_data, I_obs, UΣ)
            difs = GrIS1980s_DEM.step_through_folds(ids_train, ids_test, evaluate_fun, x_data)
            m_difs[iλ,ir]  = difs
        end
    end
    idx = vcat(ids_test...)
    # save
    to_save = (; Σ, dict, grd, λs, rs, m_difs, dists, idx=I_no_ocean[I_obs[idx]], nfiles, only_atm, method="SVD")
    dest = get_cv_file_SVD(grd, nfiles, logℓ, only_atm)
    jldsave(dest; to_save...)
end

for f in fnames_modes
    do_validation_and_save(f)
end
