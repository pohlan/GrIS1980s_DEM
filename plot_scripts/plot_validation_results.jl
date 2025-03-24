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
f_dict = [get_cv_file_SVD(grd, 70), get_cv_file_kriging(grd, maxn0)]
# plot
p_mean   = plot()
p_std    = plot()
for f in f_dict
    @unpack xc, yc, idx, grd, method = load(f)
    if method == "kriging"
        @unpack difs, binfield1, h_ref = load(f)
        dif_destd = destandardize(difs, binfield1, h_ref, add_mean=false)
        # ε = h_obs - h_rec
        # for dh:
        # dif_dh = dh_obs - dh_rec
        #        = (h_ref-h_obs) - (h_ref-h_rec) = h_rec - h_obs
        #  --> so in that case ε = - dif_dh
        dif_destd .= .- dif_destd
    elseif method == "SVD"
        @unpack λs, rs, m_difs, nfiles = load(f)
        h_ref     = h_ref_m[idx]
        iλ        = findfirst(λs .== λ0)
        ir        = findfirst(rs .== r0)
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
p_both_lin = plot(p_mean, p_std, size=(2000,700), margin=15Plots.mm, dpi=300)
savefig(joinpath(fig_dir_main, "f04.png"))


##################################################
# Kriging errors for different maxns, Figure S1  #
##################################################

# get file
fs_kriging = glob(get_cv_file_kriging(grd, "*"))
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
    @unpack idx, h_ref, difs, maxn, binfield1 = load(f)
    dif_destd = destandardize(difs, binfield1, h_ref, add_mean=false)
    dif_destd .= .- dif_destd
    dh_binned, bin_centers = svd_IceSheetDEM.bin_equal_bin_size(h_ref, dif_destd, 14)
    plot!(p_mean, bin_centers, mean.(dh_binned), label=""; color, attr...)
    hline!(p_mean, [0.0], color="grey", lw=3, z_order=1, label="", ls=:dash)
    plot!(p_std,    bin_centers, std.(dh_binned), label="$(maxn)", legend_title=L"\mathrm{n_{obs}}"; color, attr...)
end
p_mean = plot(p_mean, legend=false)
svd_IceSheetDEM.panel_annotate!(p_mean, "a")
p_std = plot(p_std, legend=:topright, legend_foreground_color=nothing)
svd_IceSheetDEM.panel_annotate!(p_std, "b")
p_both = plot(p_mean, p_std, size=(2000,700), margin=15Plots.mm, dpi=300)
savefig(joinpath(fig_dir_main, "fS01.png"))


##################################################################
# Kriging interpolated maps for different maxn values, Figure S2 #
##################################################################

Plots.scalefontsizes()
Plots.scalefontsizes(1.6)
@unpack maxns, m_interps, xsp, ysp = load(kriging_findmaxn_file())
m_interps[m_interps .== 0.] .= NaN
ps = Plots.Plot{Plots.GRBackend}[]
clims=(-4,4)
cmap = :RdBu
for m in axes(m_interps, 3)
    # cbar = m == 3 ? true : false
    yticks = m == 1 ? true : false
    pi = heatmap(x[xsp], y[ysp], m_interps[:,:,m]', title="\n"*L"\mathrm{n_{obs}} = "*"$(maxns[m])", aspect_ratio=1, wsize=(1500,400), grid=false, cbar=false, tick_direction=:out, titlefontsize=18, xlabel="Easting (m)"; clims, cmap)
    if m !== 1
        pi = plot(pi, ytickfontsize=1, ytickfontcolor=:white)
    else
        pi = plot(pi, ylabel="Northing (m)")
    end
    plot!(outl, xlims=extrema(x[xsp]), ylims=extrema(y[ysp]), fill=nothing, lw=0.5)
    push!(ps, pi)
end
p_panels = plot(ps..., size=(3000, 500), layout=(1,4), leftmargin=10Plots.mm, rightmargin=10Plots.mm, topmargin=-10Plots.mm, bottommargin=-10Plots.mm)
# colorbar
xx = range(clims...,1000)
zz = zero(xx)' .+ xx
p_c = heatmap(xx, xx, zz, ratio=15, xticks=false, legend=false, fc=cgrad(cmap), lims=clims, framestyle=:box, left_margin=-20Plots.mm, top_margin=30Plots.mm, bottom_margin=30Plots.mm, ymirror=true)
annotate!(40, 0.0, text("m", 18, "Computer Modern", :left))
# plot everything together
plot(p_panels, p_c, right_margin=[-70Plots.mm -70Plots.mm -70Plots.mm -540Plots.mm -500Plots.mm], left_margin=-20Plots.mm, bottom_margin=20Plots.mm, top_margin=10Plots.mm, size=(2400, 500), dpi=300)
savefig(joinpath(fig_dir_main, "fS02.png"))


##############################################################
# SVD error for different λ&r and singular values, Figure S3 #
##############################################################

# get file
f_SVD = f_dict[1]
# plot absolute mean of cross-validation error for different λ and r values
attr = (;size=(900,700), margin=15Plots.mm, markersize=6, lw=3.5, markerstrokewidth=0)
p_mean = plot(xlabel="Elevation of reference DEM (m)", ylabel=L"Mean $\epsilon$ (m)")
p_std  = plot(xlabel="Elevation of reference DEM (m)", ylabel=L"Std $\epsilon$ (m)")
@unpack m_difs, λs, rs, Σ = load(f_SVD)
p = plot(xscale=:log10, xlabel="λ", xticks=10 .^(5:11), ylabel="Mean absolute error (m)", title="\n Cross-validation error", size=(1300,800), grid=false,
         legend=(0.15,0.75), palette = :tol_light, legend_foreground_color=nothing)
for i in eachindex(rs)
    md_abs = [abs.(md) for md in m_difs[:,i] ]
    plot!(p, λs[1:end-1], mean.(md_abs)[1:end-1], label="r="*string(rs[i]), marker=:circle; attr...)
end
svd_IceSheetDEM.panel_annotate_xlog!(p, "a")
# plot singular values
p_sigm = plot(Σ.^2, yscale=:log10, yticks=[1e5, 1e7, 1e9, 1e11, 1e13, 1e15, 1e17], color=:black, ylabel=L"$\sigma_i^2$", xlabel=L"Mode index $i$", title="\n Singular values", label=""; attr...)
svd_IceSheetDEM.panel_annotate_ylog!(p_sigm, "b")
# save both in one figure
plot(p, p_sigm, size=(1800,600), dpi=300)
savefig(joinpath(fig_dir_main,"fS03.png"))


####################################################################
# SVD error for different number of training data files, Figure S4 #
####################################################################

# get file
fs_SVD = glob(get_cv_file_SVD(grd, "*"))
# plot
attr = (;margin=10Plots.mm, size=(svd_IceSheetDEM.wwidth,svd_IceSheetDEM.wheight), lw=1.8, markerstrokewidth=0, marker=:circle, markersize=6, markeralpha=1.0)
p_mean = plot(wsize=(svd_IceSheetDEM.wwidth, svd_IceSheetDEM.wheight))
p_std  = plot(wsize=(svd_IceSheetDEM.wwidth, svd_IceSheetDEM.wheight))
for (f_SVD, color) in zip(fs_SVD, cols)
    @unpack idx, m_difs, λs, rs, nfiles = load(f_SVD)
    if nfiles == 10   # doesn't have 500 modes
        iλ = findfirst(λs .== 1e6)
        ir = findfirst(rs .== 300)
    else
        iλ = findfirst(λs .== λ0)
        ir = findfirst(rs .== r0)
    end
    difs = m_difs[iλ, ir]
    dh_binned, bin_centers = svd_IceSheetDEM.bin_equal_bin_size(h_ref_m[idx], difs, 14)
    plot!(p_mean, bin_centers, mean.(dh_binned), legend_title=L"$m$", label=" $(nfiles*40)", ylabel=L"Mean $\epsilon$ (m)"; color, xlabel, attr...)
    hline!(p_mean, [0.0], color="grey", lw=3, z_order=1, label="", ls=:dash)
    plot!(p_std,    bin_centers, std.(dh_binned), legend_title=L"$m$", label=" $(nfiles*40)", ylabel=L"Error standard deviation $\sigma_\epsilon\quad\mathrm{(m)}$"; color, xlabel, attr...)
end
p_mean = plot(p_mean, legend_foreground_color=nothing)
svd_IceSheetDEM.panel_annotate!(p, "a")
p_std = plot(p_std, legend=false)
svd_IceSheetDEM.panel_annotate!(p_std, "b")
plot(p_mean, p_std,size=(2000,700), margin=15Plots.mm, dpi=300)
savefig(joinpath(fig_dir_main, "fS04.png"))
