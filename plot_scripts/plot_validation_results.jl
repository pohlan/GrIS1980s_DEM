using svd_IceSheetDEM, JLD2, UnPack, Glob, Plots, LaTeXStrings, Statistics, GeoStats, StatsBase, NCDatasets

# plot
Plots.scalefontsizes()
Plots.scalefontsizes(svd_IceSheetDEM.font_scaling)
attr = (;margin=10Plots.mm, size=(svd_IceSheetDEM.wwidth,svd_IceSheetDEM.wheight), lw=1.8, markerstrokewidth=0, marker=:circle, markersize=6, markeralpha=1.0)

fig_dir = joinpath("output", "main_figures")
mkpath(fig_dir)

grd = 600

f_svd = glob(joinpath("output", "validation", "cv_1e5.3_gr$(grd)_SVD_*nfiles70.jld2"))
f_krig = glob(joinpath("output", "validation", "cv_1e5.3_gr$(grd)_kriging.jld2"))
f_dict = vcat(f_svd, f_krig)

outline_shp_file = joinpath("data", "gris-imbie-1980", "gris-outline-imbie-1980_updated.shp")
csv_preprocessing, jld2_preprocessing = svd_IceSheetDEM.prepare_obs(grd, outline_shp_file)
@unpack standardize, destandardize = svd_IceSheetDEM.get_stddization_fcts(jld2_preprocessing)
@unpack href_file = load(jld2_preprocessing)

dem_ref = NCDataset(href_file)["surface"][:,:]

# choose SVD parameterization (λ and r)
λ0 = Dict("SVD_dh" => 1e7,
          "SVD_h"  => 1e7,
          "SVD_dh_detrend" => 1e5)
r0 = Dict("SVD_dh" => 500,
          "SVD_h"  => 500,
          "SVD_dh_detrend" => 500)

xlabel = "Elevation of reference DEM (m)"
cols   = Plots.palette(:tol_bright)[2:end]

p_median_all = plot()
p_std_all    = plot()
p_median_lin = plot()
p_std_lin    = plot()
p_variog = plot()

for f in f_dict
    @unpack xc, yc, idx, h_ref, grd, method = load(f)
    # if method == "SVD_h" continue end
    if method == "kriging"
        @unpack difs = load(f)
    else
        @unpack λs, rs, m_difs, nfiles = load(f)
        iλ = findfirst(λs .== λ0[method])
        ir = findfirst(rs .== r0[method])
        difs = m_difs[iλ, ir]
    end

    # calculate de-standardized difs
    if method == "SVD"
        dif_destd = difs
    else
        @unpack binfield1 = load(f)
        dif_destd = destandardize(difs, binfield1, h_ref, add_mean=false)
    end

    # ε = h_obs - h_rec
    # for dh:
    # dif_dh = dh_obs - dh_rec
    #        = (h_ref-h_obs) - (h_ref-h_rec) = h_rec - h_obs
    #  --> so in that case ε = - dif_dh
    if method == "kriging"
        dif_destd .= .- dif_destd
    end

    color = svd_IceSheetDEM.palette_dict[method]

    # heatmap
    # scatter(xc, yc, marker_z=dif_destd, cmap=:bwr, clims=(-30,30), markersize=0.7, markerstrokewidth=0, label="", aspect_ratio=1, size=(500,700))
    # savefig(joinpath(fig_dir, "eps_heatmap_"*method*".png"))

    # bin error vs elevation
    p = Plots.plot()
    dh_binned, bin_centers = svd_IceSheetDEM.bin_equal_bin_size(h_ref, dif_destd, 14)
    # boxplot
    # bc_repeat              = [zeros(length(bin_centers)).+b for b in bin_centers]
    # for (bc,dh) in zip(bc_repeat,dh_binned)
    #     boxplot!(bc, dh, label="", bar_width=220, color="cornflowerblue", markersize=1, markerstrokewidth=0, ylims=(-120,120), fillalpha=0.6, linewidth=1)
    # end
    # plot(p, ylabel=L"\epsilon_\mathrm{SVD} (m)"; xlabel, attr...)
    # hline!([0.0], color="grey", lw=3, z_order=1, label="", ls=:dot)
    # Plots.savefig(joinpath(fig_dir,"error_vs_elevation_boxplot_"*method*".png"))
    # plot only mean and std as a line plot
    # abs_dh_binned = [abs.(dh) for dh in dh_binned]
    # plot!(p_median_all, bin_centers, median.(dh_binned), label=method, ylabel=L"Median $\epsilon$ (m)"; color, xlabel, attr...)
    # hline!(p_median_all, [0.0], color="grey", lw=3, z_order=1, label="", ls=:dash)
    # plot!(p_std_all,    bin_centers, std.(dh_binned), label=method, ylabel=L"Std $\epsilon$ (m)"; color, xlabel, attr...)

    if method == "kriging" || method == "SVD"
        plot!(p_median_lin, bin_centers, mean.(dh_binned), label=split(method,"_")[1], ylabel=L"Error mean $\mu_\epsilon\quad\mathrm{(m)}$", ls=:dot; color, xlabel, attr...)
        hline!(p_median_lin, [0.0], color="grey", lw=3, z_order=1, label="", ls=:dash)
        scatter!(p_std_lin,    bin_centers, std.(dh_binned), label=split(method,"_")[1]*" bins", ylabel=L"Error standard deviation $\sigma_\epsilon\quad\mathrm{(m)}$"; color, xlabel, attr...)


        sitp, rec_errors = svd_IceSheetDEM.uncertainty_from_cv(dh_binned, bin_centers, dem_ref)
        # dest_file = joinpath(output_dir, "rec_error_$(method)_g$(grd).nc")
        # svd_IceSheetDEM.save_netcdf(dest_file, href_file, [rec_errors], ["std_error"], Dict("std_error" => Dict{String,Any}()))
        x_plot = 0.0:(diff(bin_centers)[1]*0.1):3300
        plot!(p_std_lin, x_plot, sitp.(x_plot), z_order=1, label=split(method,"_")[1]*" B-spline fit"; color)

        # variogram
        # data = svd_IceSheetDEM.make_geotable(dif_destd, xc, yc)
        # U = data |> UniqueCoords()
        # gamma = EmpiricalVariogram(U, :Z; estimator=:cressie, nlags=200,  maxlag=8e5)
        # xvals, yvals = values(gamma)
        # id_plot = findall(gamma.ordinate .!= 0.0)
        # scatter!(p_variog, xvals, yvals ./ std(dif_destd), markersize=2, label=method; attr...)
    end
end
p_median_lin = plot(p_median_lin, legend_foreground_color=nothing)
svd_IceSheetDEM.panel_annotate!(p_median_lin, "a")
p_std_lin = plot(p_std_lin, legend=:bottomleft, legend_foreground_color=nothing)
svd_IceSheetDEM.panel_annotate!(p_std_lin, "b")
p_both_lin = plot(p_median_lin, p_std_lin, size=(2000,700), margin=15Plots.mm)
savefig(joinpath(fig_dir, "error_metrics_vs_elevation.png"))

p_median_all = plot(p_median_all, legend=false)
svd_IceSheetDEM.panel_annotate!(p_median_all, "a")
p_std_all    = plot(p_std_all, legend=:bottomleft, legend_foreground_color=nothing)
svd_IceSheetDEM.panel_annotate!(p_std_all, "b")
p_both_all = plot(p_median_all, p_std_all, size=(1600,700), margin=10Plots.mm)
savefig(joinpath(fig_dir, "error_metrics_vs_elevation_all.png"))

plot(p_variog)
savefig(joinpath(fig_dir, "error_variogram.png"))


# kriging only

fs_kriging = glob(joinpath("output", "validation", "cv_1e5.3_gr$(grd)_kriging_maxn*.jld2"))

# sort files after maxn
maxns = zeros(length(fs_kriging))
for (i,f) in enumerate(fs_kriging)
    @unpack maxn = load(f)
    maxns[i]     = maxn
end
ps = sortperm(maxns)
# plot
p_median = plot(xlabel="Elevation of reference DEM (m)", ylabel=L"Median $\epsilon$ (m)")
p_std    = plot(xlabel="Elevation of reference DEM (m)", ylabel=L"Std $\epsilon$ (m)")
cols   = Plots.palette(:batlow10)[1:2:end]
for (f, color) in zip(fs_kriging[ps], cols)
    @unpack h_ref, difs, maxn, binfield1 = load(f)
    dif_destd = destandardize(difs, binfield1, h_ref, add_mean=false)
    dh_binned, bin_centers = svd_IceSheetDEM.bin_equal_bin_size(h_ref, dif_destd, 15)
    plot!(p_median, bin_centers, median.(dh_binned), label=""; color, attr...)
    hline!(p_median, [0.0], color="grey", lw=3, z_order=1, label="", ls=:dash)
    plot!(p_std,    bin_centers, std.(dh_binned), label="$(maxn)", legend_title=L"\mathrm{n_{obs}}"; color, attr...)
end
p_median = plot(p_median, legend=false)
svd_IceSheetDEM.panel_annotate!(p_median, "a")
p_std = plot(p_std, legend=:topright, legend_foreground_color=nothing)
svd_IceSheetDEM.panel_annotate!(p_std, "b")
p_both = plot(p_median, p_std, size=(2000,700), margin=15Plots.mm)
savefig(joinpath(fig_dir, "kriging_validation_maxn.png"))

# for SVD only

fs_SVD = glob(joinpath("output", "validation", "cv_*gr$(grd)*SVD_h_nfiles70_morelambdas.jld2"))

attr = (;size=(900,700), margin=15Plots.mm, markersize=6, lw=3.5, markerstrokewidth=0)

p_median = plot(xlabel="Elevation of reference DEM (m)", ylabel=L"Median $\epsilon$ (m)")
p_std    = plot(xlabel="Elevation of reference DEM (m)", ylabel=L"Std $\epsilon$ (m)")
for f_SVD in fs_SVD
    @unpack xc, yc, idx, h_ref, grd, method, m_difs, λs, rs, dict, nfiles, Σ = load(f_SVD)

    methods_name = ["mean"] #, "L2norm", "L1norm"]
    methods_fct  = [mean] #, norm, x->norm(x,1)]
    ## for SVD: plot mean, norm etc for different λ and r values
    # see https://juliagraphics.github.io/ColorSchemes.jl/stable/basics/#Pre-defined-schemes-1 for more pre-defined, color-blind friendly schemes
    for (m, fc) in zip(methods_name, methods_fct)
        p = plot(xscale=:log10, xlabel="λ", xticks=10 .^(5:11), ylabel="Mean absolute error (m)", title="\n Cross-validation error", size=(1300,800), grid=false,
                 #leftmargin=10Plots.mm, topmargin=10Plots.mm, bottommargin=10Plots.mm,
                 legend=:topleft, palette = :tol_light, legend_foreground_color=nothing)
        for i in eachindex(rs)
            md_abs = [abs.(md) for md in m_difs[:,i] ]
            plot!(p, λs[1:end-1], fc.(md_abs)[1:end-1], label="r="*string(rs[i]), marker=:circle; attr...)
        end
        svd_IceSheetDEM.panel_annotate_xlog!(p, "a")
        # plot(p)

        ########################
        # plot singular values #
        ########################
        p_sigm = plot(Σ.^2, yscale=:log10, yticks=[1e5, 1e7, 1e9, 1e11], color=:black, ylabel=L"$\sigma_i^2$", xlabel=L"Mode index $i$", title="\n Singular values", label=""; attr...)
        println(log10(ylims(p_sigm)[1]))
        println(10 .^((log10(ylims(p_sigm)[2])-log10(ylims(p_sigm)[1]))*1.07))
        println((xlims(p_sigm)[1]+(xlims(p_sigm)[2]-xlims(p_sigm)[1])*0.02, ylims(p_sigm)[1]+10 .^((log10(ylims(p_sigm)[2])-log10(ylims(p_sigm)[1]))*1.07)))
        svd_IceSheetDEM.panel_annotate_ylog!(p_sigm, "b")

        # save as subplot
        plot(p, p_sigm, size=(1800,600))
        pth = joinpath(fig_dir,m*"_abs_gr$(grd)_"*method*"_nfiles$(nfiles).png")
        println(pth)
        savefig(pth)
    end

    # plot error vs elevation for different λ and r
    # for (iλ, λ) in zip(axes(m_difs,1), λs)
    #     ir = 1
    #     dif_destd = m_difs[iλ, ir]
    #     # error vs elevation
    #     p = Plots.plot()
    #     dh_binned, bin_centers = svd_IceSheetDEM.bin_equal_bin_size(h_ref, dif_destd, 20)
    #     # plot only mean and std as a line plot
    #     md_abs = [abs.(md) for md in dh_binned ]
    #     plot!(p_median, bin_centers, median.(md_abs), label=L"$\lambda=$"*"$(λ), r=$ir"; attr...)
    #     plot!(p_std,    bin_centers, std.(dh_binned), label=L"$\lambda=$"*"$(λ), r=$ir"; attr...)
    # end
end
# plot(p_median, legend=:topright)
# plot(p_std, legend=:topright)
# savefig("p_median.png")


# plot results for different number of training files
# iλ = findfirst(λs .== 1e7)
# ir = findfirst(rs .== 250)
λ0 = 1e7
r0 = 500
fs_SVD = glob("output/validation/cv_1e5.3_gr$(grd)_SVD_h_nfiles*_morelambdas.jld2")

attr = (;margin=10Plots.mm, size=(svd_IceSheetDEM.wwidth,svd_IceSheetDEM.wheight), lw=1.8, markerstrokewidth=0, marker=:circle, markersize=6, markeralpha=1.0)

ms = []

cols = palette(:batlow10)[1:2:end]

p = plot(wsize=(svd_IceSheetDEM.wwidth, svd_IceSheetDEM.wheight))
p_std = plot(wsize=(svd_IceSheetDEM.wwidth, svd_IceSheetDEM.wheight))
for (f_SVD, color) in zip(fs_SVD, cols)
    @unpack binfield1, h_ref, m_difs, λs, rs, nfiles = load(f_SVD)
    iλ = findfirst(λs .== λ0)
    ir = findfirst(rs .== r0)
    difs = m_difs[iλ, ir]

    @unpack binfield1 = load(f_SVD)
    dif_destd = difs
    dh_binned, bin_centers = svd_IceSheetDEM.bin_equal_bin_size(h_ref, dif_destd, 14)
    # md_abs = [abs.(md) for md in dh_binned ]
    plot!(p, bin_centers, median.(dh_binned), legend_title=L"$m$", label=" $(nfiles*40)", ylabel=L"Median $\epsilon$ (m)"; color, xlabel, attr...)
    hline!(p, [0.0], color="grey", lw=3, z_order=1, label="", ls=:dash)
    plot!(p_std,    bin_centers, std.(dh_binned), legend_title=L"$m$", label=" $(nfiles*40)", ylabel=L"Error standard deviation $\sigma_\epsilon\quad\mathrm{(m)}$"; color, xlabel, attr...)
    push!(ms, norm(dif_destd,1))
end
p = plot(p, legend_foreground_color=nothing)
svd_IceSheetDEM.panel_annotate!(p, "a")
p_std = plot(p_std, legend=false)
svd_IceSheetDEM.panel_annotate!(p_std, "b")
plot(p, p_std,size=(2000,700), margin=15Plots.mm )

savefig(joinpath(fig_dir, "nfiles_comparison.png"))


##############################################################
# Figure S2: kriging interpolation for different maxn values #
##############################################################

outline_shp_file = "data/gris-imbie-1980/gris-outline-imbie-1980_updated_crs.shp"
shp              = Shapefile.shapes(Shapefile.Table(outline_shp_file))

Plots.scalefontsizes()
Plots.scalefontsizes(1.6)
@unpack maxns, m_interps = load("output/validation/kriging_findmaxn.jld2")
ps = Plots.Plot{Plots.GRBackend}[]
clims=(-4,4)
cmap = :RdBu
for m in axes(m_interps, 3)
    # cbar = m == 3 ? true : false
    yticks = m == 1 ? true : false
    pi = heatmap(x[xsp], y[ysp], m_interps[:,:,m]', title="\n"*L"\mathrm{n_{obs}} = "*"$(maxns[m])", aspect_ratio=1, wsize=(1500,400), grid=false, cbar=false, tick_direction=:out, titlefontsize=18; clims, cmap)
    if m !== 1
        pi = plot(pi, ytickfontsize=1, ytickfontcolor=:white)
    end
    plot!(shp, xlims=extrema(x[xsp]), ylims=extrema(y[ysp]), fill=nothing, lw=0.5)
    push!(ps, pi)
end
p_panels = plot(ps..., size=(3000, 500), layout=(1,4), leftmargin=10Plots.mm, rightmargin=10Plots.mm, topmargin=-10Plots.mm, bottommargin=-10Plots.mm)

xx = range(clims...,1000)
zz = zero(xx)' .+ xx
p_c = heatmap(xx, xx, zz, ratio=15, xticks=false, legend=false, fc=cgrad(cmap), lims=(-4,4), framestyle=:box, left_margin=-200Plots.mm, top_margin=30Plots.mm, bottom_margin=30Plots.mm, ymirror=true) #, size=(10,100))
annotate!(40, 0.0, text("m", 18, "Computer Modern", :left))
plot(p_c)

# plot again everything together
plot(p_panels, p_c, right_margin=-10Plots.mm) #; bottom_margin=-40Plots.mm, size=(2100,600), top_margin=10Plots.mm)
savefig("output/validation/figures/kriging_interp_maps_maxn.png")
