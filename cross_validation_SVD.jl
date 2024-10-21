using svd_IceSheetDEM, Glob, NCDatasets, JLD2, GeoStats, Statistics, StatsBase, StatsPlots, LaTeXStrings, CSV, DataFrames, UnPack, LinearAlgebra
import Plots

# set target directories
main_output_dir = joinpath("output","validation")
mkpath(main_output_dir)
fig_dir = joinpath(main_output_dir, "figures/")
mkpath(fig_dir)

# for running the script interactively
# ARGS = [
#         "--shp_file", "data/gris-imbie-1980/gris-outline-imbie-1980_updated.shp",
#         "--training_data", readdir("data/training_data_it2_600", join=true)...]

parsed_args = parse_commandline(ARGS)
model_files = parsed_args["training_data"]
outline_shp_file = parsed_args["shp_file"]
use_arpack  = parsed_args["use_arpack"]

# determine grid resolution
template_file = model_files[1]
x = NCDataset(template_file)["x"][:]
y = NCDataset(template_file)["y"][:]
const grd = Int(x[2] - x[1])

# get data
bedmachine_original, bedm_file = svd_IceSheetDEM.create_bedmachine_grid(grd)
reference_file_g150, ref_file  = svd_IceSheetDEM.create_grimpv2(grd, bedmachine_original)
h_ref        = NCDataset(ref_file)["Band1"][:]
h_ref        = svd_IceSheetDEM.replace_missing(h_ref, 0.0)

csv_preprocessing, jld2_preprocessing = svd_IceSheetDEM.prepare_obs(grd, outline_shp_file; blockspacing=grd/3, nbins1=7, nbins2=12)

# get I_no_ocean, (de-)standardization functions and variogram from pre-processing
df_all = CSV.read(csv_preprocessing, DataFrame)
dict   = load(jld2_preprocessing)
@unpack I_no_ocean, idx_aero = dict
@unpack standardize, destandardize = svd_IceSheetDEM.get_stddization_fcts(jld2_preprocessing)

# give λ and r values to loop through
λs        = [1e4, 1e5, 1e6, 1e7, 1e8]
rs        = [10, 50, 100, 150, 200, 250]

function do_validation_and_save(f)
    # load data
    @unpack U, Σ, data_mean, data_binfield1, nfiles, input = load(f)
    UΣ = U*diagm(Σ)

    # load datasets, take full SVD (to be truncated with different rs later)
    x_data, I_obs                 = svd_IceSheetDEM.prepare_obs_SVD(grd, csv_preprocessing, I_no_ocean, data_mean, main_output_dir, fig_dir; input)

    # create geotable (for GeoStats)
    x_Iobs   = x[get_ix.(I_no_ocean[I_obs],length(x))]
    y_Iobs   = y[get_iy.(I_no_ocean[I_obs],length(x))]
    geotable = svd_IceSheetDEM.make_geotable(x_data, x_Iobs, y_Iobs)

    # create sets of training and test data
    ℓ    = 2e5
    flds = folds(geotable, BlockFolding(ℓ))

    function predict_vals(λ, r, i_train, i_test, x_data, I_obs, UΣ)
        _, x_rec = svd_IceSheetDEM.solve_optim(UΣ, I_obs[i_train], r, λ, x_data[i_train])
        return x_rec[I_obs[i_test]]
    end

    # loop through λ and r values
    m_difs  = [Float32[] for i in eachindex(λs), j in eachindex(rs)]
    m_xc = [Float32[] for i in eachindex(λs), j in eachindex(rs)]
    m_yc = [Float32[] for i in eachindex(λs), j in eachindex(rs)]
    for (iλ,λ) in enumerate(λs)
        for (ir,r) in enumerate(rs)
            logλ = round(log(10, λ),digits=1)
            println("r = $r, logλ = $logλ")
            evaluate_fun(i_train,i_test) = predict_vals(λ, r, i_train, i_test, x_data, I_obs, UΣ)
            difs, xc, yc = svd_IceSheetDEM.step_through_folds(flds, evaluate_fun, geotable, save_coords=true, save_distances=false)
            m_difs[iλ,ir] = difs
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
    to_save = (; dict, grd, λs, rs, m_difs, xc=m_xc[1], yc=m_yc[1], idx, nfiles, method="SVD_"*input, binfield1=Float32.(data_binfield1[I_obs[idx]]), h_ref=Float32.(h_ref[I_no_ocean[I_obs[idx]]]))
    logℓ = round(log(10,ℓ),digits=1)
    dest = joinpath(main_output_dir,"cv_1e$(logℓ)_gr$(grd)_SVD_"*input*"_nfiles$(nfiles).jld2")
    jldsave(dest; to_save...)
end

fls = glob("output/SVD_reconstruction/SVD_components_*_nfiles*.jld2")

for f in fls
    do_validation_and_save(f)
end
