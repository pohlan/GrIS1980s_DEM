using Glob, Statistics, StatsBase

# define output path for other figures not in paper
fig_dir_others = joinpath("output", "validation", "figures")
mkpath(fig_dir_others)

# plotting attributes
Plots.scalefontsizes()
Plots.scalefontsizes(svd_IceSheetDEM.font_scaling)
attr = (;margin=10Plots.mm, size=(svd_IceSheetDEM.wwidth,svd_IceSheetDEM.wheight), lw=1.8, markerstrokewidth=0, marker=:circle, markersize=6, markeralpha=1.0)
xlabel = "Elevation of reference DEM (m)"
cols   = palette(:batlow10)[1:2:end]

##############################################
# SVD vs kriging cross-validation, Figure 4  #
##############################################

# get files
logℓ   = round(log10(2e5),digits=1)
f_dict = [get_cv_file_SVD(grd, logℓ, 70), get_cv_file_kriging(grd, logℓ, 1500)]
# choose SVD parameterization (λ and r)
λ0 = 1e7
r0 = 500
# plot
p_mean   = plot()
p_std    = plot()
for f in f_dict
    @unpack xc, yc, idx, h_ref, grd, method = load(f)
    if method == "kriging"
        @unpack difs = load(f)
        dif_destd = destandardize(difs, binfield1, h_ref, add_mean=false)
        # ε = h_obs - h_rec
        # for dh:
        # dif_dh = dh_obs - dh_rec
        #        = (h_ref-h_obs) - (h_ref-h_rec) = h_rec - h_obs
        #  --> so in that case ε = - dif_dh
        dif_destd .= .- dif_destd
    elseif method == "SVD"
        @unpack λs, rs, m_difs, nfiles = load(f)
        iλ = findfirst(λs .== λ0[method])
        ir = findfirst(rs .== r0[method])
        dif_destd = m_difs[iλ, ir]
    end
    # bin error vs elevation
    dh_binned, bin_centers = svd_IceSheetDEM.bin_equal_bin_size(h_ref, dif_destd, 14)
    # color
    color = svd_IceSheetDEM.palette_dict[method]
    # plot mean and std errors
    plot!(p_mean, bin_centers, mean.(dh_binned), label=split(method,"_")[1], ylabel=L"Error mean $\mu_\epsilon\quad\mathrm{(m)}$", ls=:dot; color, xlabel, attr...)
    hline!(p_mean, [0.0], color="grey", lw=3, z_order=1, label="", ls=:dash)
    scatter!(p_std,    bin_centers, std.(dh_binned), label=split(method,"_")[1]*" bins", ylabel=L"Error standard deviation $\sigma_\epsilon\quad\mathrm{(m)}$"; color, xlabel, attr...)
    # uncertainty B-spline fit
    sitp, rec_errors = uncertainty_from_cv(dh_binned, bin_centers, h_ref)
    x_plot = 0.0:(diff(bin_centers)[1]*0.1):3300
    plot!(p_std, x_plot, sitp.(x_plot), z_order=1, label=split(method,"_")[1]*" B-spline fit"; color)
end
p_mean = plot(p_mean, legend_foreground_color=nothing)
svd_IceSheetDEM.panel_annotate!(p_mean, "a")
p_std = plot(p_std, legend=:bottomleft, legend_foreground_color=nothing)
svd_IceSheetDEM.panel_annotate!(p_std, "b")
p_both_lin = plot(p_mean, p_std, size=(2000,700), margin=15Plots.mm)
savefig(joinpath(fig_dir_main, "Figure4.png"))


##################################################
# Kriging errors for different maxns, Figure S1  #
##################################################

# get file
fs_kriging = glob(get_cv_file_kriging(grd, logℓ, "*"))
# sort files after maxn
maxns = zeros(length(fs_kriging))
for (i,f) in enumerate(fs_kriging)
    @unpack maxn = load(f)
    maxns[i]     = maxn
end
ps = sortperm(maxns)
# plot
p_mean = plot(xlabel="Elevation of reference DEM (m)", ylabel=L"Mean $\epsilon$ (m)")
p_std  = plot(xlabel="Elevation of reference DEM (m)", ylabel=L"Std $\epsilon$ (m)")
for (f, color) in zip(fs_kriging[ps], cols)
    @unpack h_ref, difs, maxn, binfield1 = load(f)
    dif_destd = destandardize(difs, binfield1, h_ref, add_mean=false)
    dh_binned, bin_centers = svd_IceSheetDEM.bin_equal_bin_size(h_ref, dif_destd, 14)
    plot!(p_mean, bin_centers, mean.(dh_binned), label=""; color, attr...)
    hline!(p_mean, [0.0], color="grey", lw=3, z_order=1, label="", ls=:dash)
    plot!(p_std,    bin_centers, std.(dh_binned), label="$(maxn)", legend_title=L"\mathrm{n_{obs}}"; color, attr...)
end
p_mean = plot(p_mean, legend=false)
svd_IceSheetDEM.panel_annotate!(p_mean, "a")
p_std = plot(p_std, legend=:topright, legend_foreground_color=nothing)
svd_IceSheetDEM.panel_annotate!(p_std, "b")
p_both = plot(p_mean, p_std, size=(2000,700), margin=15Plots.mm)
savefig(joinpath(fig_dir_main, "FigureS1.png"))


##################################################################
# Kriging interpolated maps for different maxn values, Figure S2 #
##################################################################

Plots.scalefontsizes()
Plots.scalefontsizes(1.6)
@unpack maxns, m_interps = load("output/validation/kriging_findmaxn.jld2")   # TODO !!
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
    plot!(outl, xlims=extrema(x[xsp]), ylims=extrema(y[ysp]), fill=nothing, lw=0.5)
    push!(ps, pi)
end
p_panels = plot(ps..., size=(3000, 500), layout=(1,4), leftmargin=10Plots.mm, rightmargin=10Plots.mm, topmargin=-10Plots.mm, bottommargin=-10Plots.mm)
# colorbar
xx = range(clims...,1000)
zz = zero(xx)' .+ xx
p_c = heatmap(xx, xx, zz, ratio=15, xticks=false, legend=false, fc=cgrad(cmap), lims=(-4,4), framestyle=:box, left_margin=-200Plots.mm, top_margin=30Plots.mm, bottom_margin=30Plots.mm, ymirror=true) #, size=(10,100))
annotate!(40, 0.0, text("m", 18, "Computer Modern", :left))
plot(p_c)
# plot everything together
plot(p_panels, p_c, right_margin=-10Plots.mm)
savefig(joinpath(fig_dir_main, "FigureS2.png"))


##############################################################
# SVD error for different λ&r and singular values, Figure S3 #
##############################################################

# get file
f_SVD = f_dict[1]
# plot absolute mean of cross-validation error for different λ and r values
attr = (;size=(900,700), margin=15Plots.mm, markersize=6, lw=3.5, markerstrokewidth=0)
p_mean = plot(xlabel="Elevation of reference DEM (m)", ylabel=L"Mean $\epsilon$ (m)")
p_std  = plot(xlabel="Elevation of reference DEM (m)", ylabel=L"Std $\epsilon$ (m)")
@unpack xc, yc, idx, h_ref, grd, method, m_difs, λs, rs, dict, nfiles, Σ = load(f_SVD)
p = plot(xscale=:log10, xlabel="λ", xticks=10 .^(5:11), ylabel="Mean absolute error (m)", title="\n Cross-validation error", size=(1300,800), grid=false,
         legend=:topleft, palette = :tol_light, legend_foreground_color=nothing)
for i in eachindex(rs)
    md_abs = [abs.(md) for md in m_difs[:,i] ]
    plot!(p, λs[1:end-1], mean.(md_abs)[1:end-1], label="r="*string(rs[i]), marker=:circle; attr...)
end
svd_IceSheetDEM.panel_annotate_xlog!(p, "a")
# plot singular values
p_sigm = plot(Σ.^2, yscale=:log10, yticks=[1e5, 1e7, 1e9, 1e11], color=:black, ylabel=L"$\sigma_i^2$", xlabel=L"Mode index $i$", title="\n Singular values", label=""; attr...)
svd_IceSheetDEM.panel_annotate_ylog!(p_sigm, "b")
# save both in one figure
plot(p, p_sigm, size=(1800,600))
savefig(joinpath(fig_dir_main,"FigureS3.png"))


####################################################################
# SVD error for different number of training data files, Figure S4 #
####################################################################

# get file
fs_SVD = glob(get_cv_file_SVD(grd, logℓ, "*"))
# plot
attr = (;margin=10Plots.mm, size=(svd_IceSheetDEM.wwidth,svd_IceSheetDEM.wheight), lw=1.8, markerstrokewidth=0, marker=:circle, markersize=6, markeralpha=1.0)
p_mean = plot(wsize=(svd_IceSheetDEM.wwidth, svd_IceSheetDEM.wheight))
p_std  = plot(wsize=(svd_IceSheetDEM.wwidth, svd_IceSheetDEM.wheight))
for (f_SVD, color) in zip(fs_SVD, cols)
    @unpack binfield1, h_ref, m_difs, λs, rs, nfiles = load(f_SVD)
    iλ = findfirst(λs .== λ0)
    ir = findfirst(rs .== r0)
    difs = m_difs[iλ, ir]
    dh_binned, bin_centers = svd_IceSheetDEM.bin_equal_bin_size(h_ref, difs, 14)
    plot!(p_mean, bin_centers, mean.(dh_binned), legend_title=L"$m$", label=" $(nfiles*40)", ylabel=L"Mean $\epsilon$ (m)"; color, xlabel, attr...)
    hline!(p_mean, [0.0], color="grey", lw=3, z_order=1, label="", ls=:dash)
    plot!(p_std,    bin_centers, std.(dh_binned), legend_title=L"$m$", label=" $(nfiles*40)", ylabel=L"Error standard deviation $\sigma_\epsilon\quad\mathrm{(m)}$"; color, xlabel, attr...)
end
p_mean = plot(p_mean, legend_foreground_color=nothing)
svd_IceSheetDEM.panel_annotate!(p, "a")
p_std = plot(p_std, legend=false)
svd_IceSheetDEM.panel_annotate!(p_std, "b")
plot(p_mean, p_std,size=(2000,700), margin=15Plots.mm )
savefig(joinpath(fig_dir_main, "FigureS4.png"))
