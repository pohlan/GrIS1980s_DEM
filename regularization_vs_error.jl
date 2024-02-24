using svd_IceSheetDEM, Glob, NCDatasets, JLD2
import Plots

use_arpack = false

shp_file = "data/gris-imbie-1980/gris-outline-imbie-1980_updated.shp"
model_files   = glob("data/training_data_it2_600/usurf*.nc")[1:2]
template_file = model_files[1]
x = NCDataset(template_file)["x"]
y = NCDataset(template_file)["y"]
const gr = Int(x[2] - x[1])
const F  = Float32

bedm_file              = create_bedmachine_grid(gr, template_file)
bedmachine_path        = splitdir(bedm_file)[1]
aerodem_g150, obs_file = create_aerodem(;gr, shp_file, bedmachine_path)
imbie_mask             = create_imbie_mask(;gr, shp_file, sample_path=aerodem_g150)

fig_dir = "output/figures/"
mkpath(fig_dir)
destdir = "output/model_selection/"
mkpath(destdir)

r_gap  = 2e4
n_gaps = 2

λs = [1e5, 1e6, 1e7, 5e7, 1e8, 5e8, 1e9, 1e10]
rs = [60, 80, 100, 140, 160, 200, 250, 300, 350]

# load datasets, take full SVD (to be truncated with different rs later)
UΣ, I_no_ocean, Data_mean, Σ = svd_IceSheetDEM.prepare_model(model_files, imbie_mask, bedm_file, F, use_arpack)
x_data, I_obs                = svd_IceSheetDEM.prepare_obs(obs_file, I_no_ocean, Data_mean)
f_eval(λ, r, i_train       ) = svd_IceSheetDEM.solve_optim(UΣ, I_obs[i_train], r, λ, x_data[i_train])

# create gaps
i_train_sets, i_test_sets    = svd_IceSheetDEM.make_training_gaps(x, y, r_gap, n_gaps, I_no_ocean, I_obs)
dict                         = svd_IceSheetDEM.sample_param_space(f_eval, λs, rs, x_data, I_obs, i_train_sets, i_test_sets)
dict["σ"] = Σ
dest = destdir*"dict_rgap_$(r_gap)_ngaps_$(n_gaps).jld2"
save(dest, dict)

# plot
for m in keys(filter(p -> .!(p.first ∈ ["r","λ", "σ"]), dict))
    Plots.plot(λs, abs.(dict[m][:,1]), xscale=:log10, xlabel="λ", ylabel=m, label="r="*string(rs[1]), marker=:circle, markerstrokewidth=0)
    for i in 2:length(rs)
        Plots.plot!(λs, abs.(dict[m][:,i]), xscale=:log10, xlabel="λ", ylabel=m, label="r="*string(rs[i]), marker=:circle, markerstrokewidth=0)
    end
    Plots.savefig(fig_dir*m*"_abs_lineplot_rgap_$(r_gap)_.png")
end
