using svd_IceSheetDEM, Glob, NCDatasets, JLD2, GeoStats, Statistics, StatsBase
import Plots

# for running the script interactively
# ARGS = [
#         "--lambda", "1e5",
#         "--r", "377",
#         "--shp_file", "data/gris-imbie-1980/gris-outline-imbie-1980_updated.shp",
#         "--training_data", readdir("data/training_data_it2_600", join=true)...]

get_ix(i,nx) = i % nx == 0 ? nx : i % nx
get_iy(i,nx) = cld(i,nx)

parsed_args = parse_commandline(ARGS)
model_files = parsed_args["training_data"]
shp_file    = parsed_args["shp_file"]
use_arpack  = parsed_args["use_arpack"]

template_file = model_files[1]
x = NCDataset(template_file)["x"][:]
y = NCDataset(template_file)["y"][:]
const gr = Int(x[2] - x[1])
const F  = Float32

bedm_file              = create_bedmachine_grid(gr, template_file)
bedmachine_path        = splitdir(bedm_file)[1]
aerodem_g150, obs_file = create_aerodem(;gr, shp_file, bedmachine_path)
imbie_mask             = create_imbie_mask(;gr, shp_file, sample_path=aerodem_g150)

fig_dir = "output/model_selection/figures/"
mkpath(fig_dir)
destdir = "output/model_selection/"
mkpath(destdir)

# give λ and r values to loop through
λs        = [1e4, 1e5] #, 1e6, 5e6, 1e7, 5e7, 1e8, 5e8, 1e9, 5e9]
rs        = [50, 100] #, 150, 200, 250, 300, 350]

# load datasets, take full SVD (to be truncated with different rs later)
UΣ, I_no_ocean, Data_mean, Σ = svd_IceSheetDEM.prepare_model(model_files[1:3], imbie_mask, bedm_file, maximum(rs), F, use_arpack)
x_data, I_obs                = svd_IceSheetDEM.prepare_obs(obs_file, I_no_ocean, Data_mean)

function evaluate_params(λ, r, i_train, i_test, x_data, I_obs, UΣ)
    x_rec = svd_IceSheetDEM.solve_optim(UΣ, I_obs[i_train], r, λ, x_data[i_train])
    return x_rec[I_obs[i_test]]
end

# create geotable (for GeoStats)
x_Iobs   = x[get_ix.(I_no_ocean[I_obs],length(x))]
y_Iobs   = y[get_iy.(I_no_ocean[I_obs],length(x))]
table    = (; Z=Float32.(x_data))
coords   = [(xi,yi) for (xi,yi) in zip(x_Iobs, y_Iobs) ]
geotable = georef(table,coords)


# create sets of training and test data
ℓ    = 7e5
flds = folds(geotable, BlockFolding(ℓ))

# loop through λ and r values
methods_name = ["median", "mean", "nmad", "std", "L2norm"]
methods_fct  = [median, mean, mad, std, norm]
# initiate dictionary
dict = Dict{String,Array}(n => zeros(length(λs), length(rs)) for n in methods_name)
dict["λ"] = λs
dict["r"] = rs
for (iλ,λ) in enumerate(λs)
    for (ir,r) in enumerate(rs)
        logλ = round(log(10, λ),digits=1)
        println("r = $r, logλ = $logλ")
        evaluate_fun(i_train,i_test) = evaluate_params(λ, r, i_train, i_test, x_data, I_obs, UΣ)
        difs = svd_IceSheetDEM.step_through_folds(flds, evaluate_fun,geotable)
        for (mn, mf) in zip(methods_name, methods_fct)
            dict[mn][iλ,ir] = mf(difs)
        end
    end
end
dict["σ"] = Σ
dest      = destdir*"dict_cv_block_$(ℓ).jld2"
save(dest, dict)

# plot
Plots.scalefontsizes()
Plots.scalefontsizes(2.5)
# see https://juliagraphics.github.io/ColorSchemes.jl/stable/basics/#Pre-defined-schemes-1 for more pre-defined, color-blind friendly schemes
for m in keys(filter(p -> .!(p.first ∈ ["r","λ", "σ"]), dict))
    p = Plots.plot(xscale=:log10, xlabel="λ", ylabel=m, size=(1300,800), leftmargin=10Plots.mm, topmargin=10Plots.mm, bottommargin=10Plots.mm, legend=:top, palette = :tol_light)
    for i in eachindex(rs)
        Plots.plot!(λs, abs.(dict[m][:,i]), label="r="*string(rs[i]), marker=:circle, markersize=6, markerstrokewidth=0, lw=3.5)
    end
    Plots.plot(p)
    logℓ = round(log(10,ℓ),digits=1)
    Plots.savefig(fig_dir*m*"_abs_cv_block_1e$(logℓ).png")
end
