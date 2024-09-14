using svd_IceSheetDEM, Glob, NCDatasets, JLD2, GeoStats, Statistics, StatsBase, StatsPlots, LaTeXStrings, CSV, DataFrames, UnPack
import Plots

# set target directories
main_output_dir = joinpath("output","validation")
mkpath(main_output_dir)
fig_dir = joinpath(main_output_dir, "figures/")
mkpath(fig_dir)

# for running the script interactively
ARGS = [
        "--shp_file", "data/gris-imbie-1980/gris-outline-imbie-1980_updated.shp",
        "--training_data", readdir("data/training_data_it1_1200", join=true)...]

parsed_args = parse_commandline(ARGS)
model_files = parsed_args["training_data"]
outline_shp_file = parsed_args["shp_file"]
use_arpack  = parsed_args["use_arpack"]

# determine grid resolution
template_file = model_files[1]
x = NCDataset(template_file)["x"][:]
y = NCDataset(template_file)["y"][:]
const gr = Int(x[2] - x[1])
const F  = Float32

# get data
bedmachine_original, bedm_file = svd_IceSheetDEM.create_bedmachine_grid(gr)
reference_file_g150, ref_file  = svd_IceSheetDEM.create_grimpv2(gr, bedmachine_original)
h_ref        = NCDataset(ref_file)["Band1"][:]
h_ref        = svd_IceSheetDEM.replace_missing(h_ref, 0.0)

csv_preprocessing, jld2_preprocessing = svd_IceSheetDEM.prepare_obs(gr, outline_shp_file; blockspacing=gr/3, nbins1=7, nbins2=12)

# get I_no_ocean, (de-)standardization functions and variogram from pre-processing
df_all = CSV.read(csv_preprocessing, DataFrame)
dict   = load(jld2_preprocessing)
@unpack I_no_ocean, idx_aero, params = dict
@unpack standardize, destandardize = svd_IceSheetDEM.get_stddization_fcts(jld2_preprocessing)

# give λ and r values to loop through
λs        = [5e6, 1e7, 5e7, 1e8]
rs        = [150, 200, 250]

# load datasets, take full SVD (to be truncated with different rs later)
x_data, I_obs                 = svd_IceSheetDEM.prepare_obs(gr, csv_preprocessing, I_no_ocean, fig_dir)
UΣ, data_mean, data_ref, Σ, _ = svd_IceSheetDEM.prepare_model(model_files, standardize, h_ref, I_no_ocean, maximum(rs), use_arpack, main_output_dir) # read in model data and take svd to derive "eigen ice sheets"

function predict_vals(λ, r, i_train, i_test, x_data, I_obs, UΣ)
    _, x_rec = svd_IceSheetDEM.solve_optim(UΣ, I_obs[i_train], r, λ, x_data[i_train])
    return x_rec[I_obs[i_test]]
end

# create geotable (for GeoStats)
x_Iobs   = x[get_ix.(I_no_ocean[I_obs],length(x))]
y_Iobs   = y[get_iy.(I_no_ocean[I_obs],length(x))]
geotable = svd_IceSheetDEM.make_geotable(x_data, x_Iobs, y_Iobs)

# create sets of training and test data
ℓ    = 2e5
flds = folds(geotable, BlockFolding(ℓ))

# loop through λ and r values
methods_name = ["median", "mean", "nmad", "std", "L2norm"]
methods_fct  = [median, mean, mad, std, norm]
dict = Dict{String,Any}(n => zeros(length(λs), length(rs)) for n in methods_name)
m_difs  = [Float32[] for i in eachindex(λs), j in eachindex(rs)]
# m_dists = [Float32[] for i in eachindex(λs), j in eachindex(rs)]
m_xc = [Float32[] for i in eachindex(λs), j in eachindex(rs)]
m_yc = [Float32[] for i in eachindex(λs), j in eachindex(rs)]
for (iλ,λ) in enumerate(λs)
    for (ir,r) in enumerate(rs)
        logλ = round(log(10, λ),digits=1)
        println("r = $r, logλ = $logλ")
        evaluate_fun(i_train,i_test) = predict_vals(λ, r, i_train, i_test, x_data, I_obs, UΣ)
        difs, xc, yc = svd_IceSheetDEM.step_through_folds(flds, evaluate_fun, geotable, save_coords=true, save_distances=false)
        # , dists, xc, yc
        for (mn, mf) in zip(methods_name, methods_fct)
            dict[mn][iλ,ir] = mf(difs)
        end
        m_difs[iλ,ir] = difs
        # m_dists[iλ,ir] = dists
        m_xc[iλ,ir] = xc
        m_yc[iλ,ir] = yc
    end
end

idxs = [Int[] for i in flds]
for (i,fs) in enumerate(flds)
    idxs[i] = fs[2]
end
idx = vcat(idxs...)

# save
to_save = (; dict, gr, λs, rs, m_difs, xc=m_xc[1], yc=m_yc[1], idx, method="SVD", h_ref=Float32.(h_ref[I_no_ocean[I_obs[idx]]]))
logℓ = round(log(10,ℓ),digits=1)
dest = joinpath(main_output_dir,"dict_cv_block_1e$(logℓ)_gr$(gr)_SVD.jld2")
jldsave(dest; to_save...)
