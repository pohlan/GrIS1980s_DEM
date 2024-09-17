using svd_IceSheetDEM, JLD2, UnPack, Glob, Plots, LaTeXStrings, Statistics, GeoStats, StatsBase

# plot
Plots.scalefontsizes()
Plots.scalefontsizes(1.9)
attr = (;margin=10Plots.mm, size=(1000,700), lw=2, markerstrokewidth=0, marker=:circle)

output_dir = "output/validation/"
fig_dir    = joinpath(output_dir, "figures/")
mkpath(fig_dir)

grd = 1200

f_dict = glob("output/validation/dict_cv_block_1e5.3_gr$(grd)_*.jld2")

outline_shp_file = "data/gris-imbie-1980/gris-outline-imbie-1980_updated.shp"
csv_preprocessing, jld2_preprocessing = svd_IceSheetDEM.prepare_obs(grd, outline_shp_file)
@unpack standardize, destandardize = svd_IceSheetDEM.get_stddization_fcts(jld2_preprocessing)

# choose SVD parameterization (λ and r)
i_SVD = 3

p_median = plot()
p_std    = plot()
p_variog = plot()
for f in f_dict
    @unpack xc, yc, idx, h_ref, gr, method = load(f)
    # if method == "SVD_h" continue end
    if method == "kriging"
        @unpack difs = load(f)
    else
        difs = m_difs[i_SVD]
    end

    # calculate de-standardized difs
    if method == "SVD_non-standardized" || method =="SVD_h"
        dif_destd = difs
    else
        dif_destd = destandardize(difs, h_ref, add_mean=false)
    end

    # heatmap
    scatter(xc, yc, marker_z=dif_destd, cmap=:bwr, clims=(-4,4), markersize=0.7, markerstrokewidth=0, label="", aspect_ratio=1, size=(500,700))
    savefig(joinpath(fig_dir, "eps_heatmap_"*method*".png"))

    # error vs elevation
    p = Plots.plot()
    dh_binned, bin_centers = svd_IceSheetDEM.bin_equal_bin_size(h_ref, dif_destd, 20)
    bc_repeat              = [zeros(length(bin_centers)).+b for b in bin_centers]
    for (bc,dh) in zip(bc_repeat,dh_binned)
        boxplot!(bc, dh, label="", bar_width=220, color="cornflowerblue", markersize=1, markerstrokewidth=0, ylims=(-120,120), fillalpha=0.6, linewidth=1)
    end
    plot(p, xlabel="Elevation of reference DEM (m)", ylabel=L"\epsilon_\mathrm{SVD} (m)"; attr...)
    hline!([0.0], color="grey", lw=3, z_order=1, label="", ls=:dot)
    Plots.savefig(joinpath(fig_dir,"error_vs_elevation_boxplot_"*method*".png"))
    # plot only mean and std as a line plot
    plot!(p_median, bin_centers, abs.(median.(dh_binned)), label=method, xlabel="Elevation of reference DEM (m)", ylabel=L"Median $\epsilon$ (m)", yscale=:log10; attr...)
    plot!(p_std,    bin_centers, std.(dh_binned), label=method, xlabel="Elevation of reference DEM (m)", ylabel=L"Std $\epsilon$ (m)", yscale=:log10; attr...)

    # variogram
    data = svd_IceSheetDEM.make_geotable(dif_destd, xc, yc)
    U = data |> UniqueCoords()
    gamma = EmpiricalVariogram(U, :Z; estimator=:cressie, nlags=1000,  maxlag=1e6)
    id_plot = findall(gamma.ordinate .!= 0.0)
    scatter!(p_variog, gamma.abscissa[id_plot], gamma.ordinate[id_plot], markersize=2, label=method, yscale=:log10; attr...)
end
p_median = plot(p_median, legend=false)
p_std = plot(p_std, legend=:bottomleft, legend_foreground_color=nothing)
p_both = plot(p_median, p_std, size=(1600,700), margin=10Plots.mm)
savefig(joinpath(fig_dir, "error_metrics_vs_elevation.png"))
plot(p_variog)
savefig(joinpath(fig_dir, "error_variogram.png"))


# for SVD only

f_SVD = "output/validation/*gr$(grd)_*_SVD_dh_detrend.jld2"

@unpack xc, yc, idx, h_ref, gr, method, m_difs, λs, rs, dict = load(f_SVD)

methods_name = ["median", "mean", "nmad", "std", "L2norm"]
methods_fct  = [median, mean, mad, std, norm]
## for SVD: plot mean, norm etc for different λ and r values
# see https://juliagraphics.github.io/ColorSchemes.jl/stable/basics/#Pre-defined-schemes-1 for more pre-defined, color-blind friendly schemes
for (m, fc) in zip(methods_name, methods_fct)
    p = plot(xscale=:log10, xlabel="λ", ylabel=m, size=(1300,800), leftmargin=10Plots.mm, topmargin=10Plots.mm, bottommargin=10Plots.mm, legend=:top, palette = :tol_light)
    for i in eachindex(rs)
        plot!(λs, abs.(fc.(m_difs[:,i])), label="r="*string(rs[i]), marker=:circle, markersize=6, markerstrokewidth=0, lw=3.5)
    end
    plot(p)
    savefig(joinpath(fig_dir,m*"_abs_gr$(gr)_"*method*".png"))
end
