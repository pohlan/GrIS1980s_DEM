using svd_IceSheetDEM, JLD2, UnPack, Glob, Plots, LaTeXStrings, Statistics, GeoStats, StatsBase

# plot
Plots.scalefontsizes()
Plots.scalefontsizes(1.9)
attr = (;margin=10Plots.mm, size=(1000,700), lw=2, markerstrokewidth=0, marker=:circle)

output_dir = "output/validation/"
fig_dir    = joinpath(output_dir, "figures/")
mkpath(fig_dir)

grd = 600

f_svd = glob("output/validation/cv_1e5.3_gr$(grd)_SVD_*nfiles70.jld2")
f_krig = glob("output/validation/cv_1e5.3_gr$(grd)_kriging.jld2")
f_dict = vcat(f_svd, f_krig)

outline_shp_file = "data/gris-imbie-1980/gris-outline-imbie-1980_updated.shp"
csv_preprocessing, jld2_preprocessing = svd_IceSheetDEM.prepare_obs(grd, outline_shp_file)
@unpack standardize, destandardize = svd_IceSheetDEM.get_stddization_fcts(jld2_preprocessing)
dict = load(jld2_preprocessing)
@unpack I_no_ocean = dict

# choose SVD parameterization (λ and r)
λ0 = 1e7
r0 = 300

xlabel = "Elevation of reference DEM (m)"
cols   = Plots.palette(:tol_bright)[2:end]

p_median_all = plot()
p_std_all    = plot()
p_median_lin = plot()
p_std_lin    = plot()
p_variog = plot()

for (f,color) in zip(f_dict, cols)
    @unpack xc, yc, idx, h_ref, grd, method = load(f)
    # if method == "SVD_h" continue end
    if method == "kriging"
        @unpack difs = load(f)
    else
        @unpack λs, rs, m_difs, nfiles = load(f)
        iλ = findfirst(λs .== λ0)
        ir = findfirst(rs .== r0)
        difs = m_difs[iλ, ir]
    end

    # calculate de-standardized difs
    if method == "SVD_dh" || method =="SVD_h"
        dif_destd = difs
    else
        @unpack binfield1 = load(f)
        dif_destd = destandardize(difs, binfield1, h_ref, add_mean=false)
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
    plot(p, ylabel=L"\epsilon_\mathrm{SVD} (m)"; xlabel, attr...)
    hline!([0.0], color="grey", lw=3, z_order=1, label="", ls=:dot)
    Plots.savefig(joinpath(fig_dir,"error_vs_elevation_boxplot_"*method*".png"))
    # plot only mean and std as a line plot
    abs_dh_binned = [abs.(dh) for dh in dh_binned]
    plot!(p_median_all, bin_centers, median.(abs_dh_binned), label=method, ylabel=L"Median $|\epsilon|$ (m)", yscale=:log10; color, xlabel, attr...)
    plot!(p_std_all,    bin_centers, std.(dh_binned), label=method, ylabel=L"Std $\epsilon$ (m)", yscale=:log10; color, xlabel, attr...)

    if method == "kriging" || method == "SVD_dh_detrend"
        plot!(p_median_lin, bin_centers, median.(dh_binned), label=method, ylabel=L"Median $\epsilon$ (m)"; color, xlabel, attr...)
        hline!(p_median_lin, [0.0], color="grey", lw=3, z_order=1, label="", ls=:dash)
        plot!(p_std_lin,    bin_centers, std.(dh_binned), label=method, ylabel=L"Std $\epsilon$ (m)"; color, xlabel, attr...)

        # variogram
        data = svd_IceSheetDEM.make_geotable(dif_destd, xc, yc)
        U = data |> UniqueCoords()
        gamma = EmpiricalVariogram(U, :Z; estimator=:cressie, nlags=200,  maxlag=5e5)
        id_plot = findall(gamma.ordinate .!= 0.0)
        scatter!(p_variog, gamma.abscissa[id_plot], gamma.ordinate[id_plot], markersize=2, label=method; attr...)
    end
end
p_median_all = plot(p_median_all, legend=false)
p_std_all    = plot(p_std_all, legend=:bottomleft, legend_foreground_color=nothing)
p_both_all = plot(p_median_all, p_std_all, size=(1600,700), margin=10Plots.mm)
savefig(joinpath(fig_dir, "error_metrics_vs_elevation_all_logscale.png"))
p_median_lin = plot(p_median_lin, legend=false)
p_std_lin = plot(p_std_lin, legend=:bottomleft, legend_foreground_color=nothing)
p_both_lin = plot(p_median_lin, p_std_lin, size=(1600,700), margin=10Plots.mm)
savefig(joinpath(fig_dir, "error_metrics_vs_elevation.png"))
plot(p_variog)
savefig(joinpath(fig_dir, "error_variogram.png"))


# for SVD only

fs_SVD = glob("output/validation/cv_*gr$(grd)*SVD_dh_detrend_nfiles70.jld2")

for f_SVD in fs_SVD
    @unpack xc, yc, idx, h_ref, grd, method, m_difs, λs, rs, dict, nfiles = load(f_SVD)

    methods_name = ["mean"] #, "L2norm", "L1norm"]
    methods_fct  = [mean] #, norm, x->norm(x,1)]
    ## for SVD: plot mean, norm etc for different λ and r values
    # see https://juliagraphics.github.io/ColorSchemes.jl/stable/basics/#Pre-defined-schemes-1 for more pre-defined, color-blind friendly schemes
    for (m, fc) in zip(methods_name, methods_fct)
        p = plot(xscale=:log10, xlabel="λ", ylabel="mean absolute error (m)", size=(1300,800), leftmargin=10Plots.mm, topmargin=10Plots.mm, bottommargin=10Plots.mm, legend=:top, palette = :tol_light)
        for i in eachindex(rs)
            md_abs = [abs.(md) for md in m_difs[:,i] ]
            plot!(λs, fc.(md_abs), label="r="*string(rs[i]), marker=:circle, markersize=6, markerstrokewidth=0, lw=3.5)
        end
        plot(p)
        savefig(joinpath(fig_dir,m*"_abs_gr$(grd)_"*method*"_nfiles$(nfiles).png"))
    end
end

# plot error vs elevation for different λ and r
# p_median = plot(xlabel="Elevation of reference DEM (m)", ylabel=L"Median $\epsilon$ (m)")
# p_std    = plot(xlabel="Elevation of reference DEM (m)", ylabel=L"Std $\epsilon$ (m)")
# for (iλ, λ) in zip(axes(m_difs,1), λs)
#     # for (ir, r) in zip(axes(m_difs,2), rs)
#     ir = 3
#     dif_destd = m_difs[iλ, ir]
#     # error vs elevation
#     p = Plots.plot()
#     dh_binned, bin_centers = svd_IceSheetDEM.bin_equal_bin_size(h_ref, dif_destd, 20)
#     # plot only mean and std as a line plot
#     md_abs = [abs.(md) for md in dh_binned ]
#     plot!(p_median, bin_centers, norm.(md_abs), label=L"$\lambda=$"*"$(λ), r=$ir"; attr...)
#     plot!(p_std,    bin_centers, std.(dh_binned), label=L"$\lambda=$"*"$(λ), r=$ir"; attr...)
#     # end
# end
# plot(p_median, legend=:topright)



# plot results for different number of training files
# iλ = findfirst(λs .== 1e7)
# ir = findfirst(rs .== 250)
iλ = 1
ir = 1
fs_SVD = glob("output/validation/cv_1e5.3_gr$(grd)_SVD_dh_detrend_*nfiles*.jld2")

ms = []
p = plot()
for f_SVD in fs_SVD[2:end]
    @unpack binfield1, h_ref, m_difs, λs, rs, nfiles = load(f_SVD)
    difs = m_difs[iλ, ir]

    @unpack binfield1 = load(f_SVD)
    dif_destd = destandardize(difs, binfield1, h_ref, add_mean=false)

    dh_binned, bin_centers = svd_IceSheetDEM.bin_equal_bin_size(h_ref, dif_destd, 20)
    md_abs = [abs.(md) for md in dh_binned ]
    plot!(p, bin_centers, median.(dh_binned), label="nfiles=$(nfiles)", ylabel=L"Median $\epsilon$ (m)"; xlabel, attr...)
    hline!(p, [0.0], color="grey", lw=3, z_order=1, label="", ls=:dash)
    push!(ms, norm(dif_destd,1))
end
plot(p)
savefig(joinpath(fig_dir, "nfiles_comparison.png"))

