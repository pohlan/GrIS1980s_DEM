# define output path for other figures not in paper
fig_dir_others = joinpath("output", "validation", "figures")
mkpath(fig_dir_others)

# plotting attributes
Plots.scalefontsizes()
Plots.scalefontsizes(GrIS1980s_DEM.font_scaling)
attr = (;margin=10Plots.mm, size=(GrIS1980s_DEM.wwidth,GrIS1980s_DEM.wheight), lw=1.8, markerstrokewidth=0, marker=:circle, markersize=6, markeralpha=1.0)

##############################################
# SVD vs kriging cross-validation, Figure 4  #
##############################################

function make_validation_plots(l1, l2, l3; standardized::Bool)
    f_dict = [get_cv_file_GP(grd, only_atm=true), get_cv_file_SVD(grd, 70, only_atm=false)]
    f_SVD_only_atm = get_cv_file_SVD(grd, 70, only_atm=true)
    p_mean   = plot()
    p_std    = plot()
    p_box    = plot()
    for f in f_dict
        @unpack idx, grd, method = load(f)
        if method == "GP"
            @unpack difs, dists, sigmas, binfield1, h_ref = load(f)
            # ε = h_obs - h_rec
            # for dh:
            # dif_dh = dh_obs - dh_rec
            #        = (h_ref-h_obs) - (h_ref-h_rec) = h_rec - h_obs
            #  --> so in that case ε = - dif_dh
            difs .= .- difs
            if !standardized
                difs = destandardize(difs, binfield1, h_ref, add_mean=false)
                bin_metric = h_ref
            else
                bin_metric = vcat(dists...)
            end
        elseif method == "SVD"
            @unpack λs, rs, m_difs, dists, nfiles = load(f)
            iλ         = findfirst(λs .== λ0)
            ir         = findfirst(rs .== r0)
            difs       = m_difs[iλ, ir]
            if standardized
                # difs       = standardize(difs, bfield_1_m[idx], h_ref_m[idx], subtract_mean=false)
                bin_metric = vcat(dists...)
                ii = findall(bin_metric .> 110e3)
                deleteat!(bin_metric, ii)
                deleteat!(difs, ii)
            else
                bin_metric = h_ref_m[idx]
            end
        end
        # bin error vs elevation
        dh_binned, bin_centers = GrIS1980s_DEM.bin_equal_bin_size(bin_metric, difs, 14)
        # color
        color = GrIS1980s_DEM.palette_dict[method]
        # plot mean and std errors
        label = method == "GP" ? split(method,"_")[1]*" (only ATM)" : split(method,"_")[1]*"(all)"
        ylabel_mean = standardized ? L"Error mean $\mu_\epsilon$ (standardized)" : L"Error mean $\mu_\epsilon\quad\mathrm{(m)}$"
        xlabel = standardized ? "Distance to closest observation (km)" : "Elevation of reference DEM (m)"
        xscl    = standardized ? 1e-3 : 1
        if !standardized
            plot!(p_std, bin_centers.*xscl, std.(dh_binned), ylabel=L"Error STD $\sigma_\epsilon\quad\mathrm{(m)}$", ls=:dot, ylims=(0,43); color, xlabel, label, attr...)
            bar!(p_std, bin_centers.*xscl, length.(dh_binned)./1600, alpha=0.4; color, label=method*" rel. bin size")
            plot!(p_mean, bin_centers.*xscl, mean.(dh_binned), ylabel=ylabel_mean, ls=:dot; color, xlabel, label, attr...)
            hline!(p_mean, [0.0], color="grey", lw=3, z_order=1, label="", ls=:dash) #, alpha=0.7; color)
            if method =="SVD" GrIS1980s_DEM.panel_annotate!(p_mean, l2) end
        elseif method == "GP"
            plot!(p_std, bin_centers.*xscl, std.(dh_binned), ylabel=L"Error STD $\sigma_\epsilon$ (standardized)", ls=:dot, y_foreground_color_text=color, y_guidefontcolor=color; label="", color, xlabel, attr...)
            plot!(p_mean, bin_centers.*xscl, mean.(dh_binned), ylabel=ylabel_mean, ls=:dot, y_foreground_color_text=color, y_guidefontcolor=color, legend=false; label="", color, xlabel, attr...)
            hline!(p_mean, [0.0], lw=3, z_order=1, label="", ls=:dash, alpha=0.7; color)
            GrIS1980s_DEM.panel_annotate!(p_mean, l2)
            sig_binned, sig_bins = GrIS1980s_DEM.bin_equal_bin_size(bin_metric, sqrt.(sigmas), 14)
            plot!(p_std, sig_bins.*xscl, mean.(sig_binned).- 0.35, color=:cyan2, marker=:x, z_order=1, label="GP uncertainty") #, alpha=0.4)
            # plot!(p_std, sig_bins.*xscl, max.(mean.(sig_binned).- 2 .*std.(sig_binned).- 0.35, 0.0), fillrange=mean.(sig_binned).+ 2 .*std.(sig_binned).-0.35, color=:deepskyblue2, alpha=0.4)
            # histogram2d!(p_std, bin_metric.*xscl, sigmas .- std(dh_binned[1]), bins=180, cmap=:dense, colorbar=false)
            bar!(p_std, bin_centers.*xscl, length.(dh_binned)./75000, alpha=0.3, label=""; color)
         else
            p_tw_m = twinx(p_mean)
            plot!(p_tw_m, bin_centers.*xscl, mean.(dh_binned), ylabel=L"Error mean $\mu_\epsilon\quad\mathrm{(m)}$", ls=:dot, y_foreground_color_text=color, y_guidefontcolor=color, legend=false; color, attr...)
            hline!(p_tw_m, [0.0], lw=3, z_order=1, label="", ls=:dash, alpha=0.7; color)
            p_tw_std = twinx(p_std)
            plot!(p_tw_std, bin_centers.*xscl, std.(dh_binned), ylabel=L"Error STD $\sigma_\epsilon\quad\mathrm{(m)}$", ls=:dot, y_foreground_color_text=color, y_guidefontcolor=color, label=""; color, attr...)
            bar!(p_std, bin_centers.*xscl, length.(dh_binned)./75000, alpha=0.3, label=""; color)
        end
        boxplot!(p_box, difs, outliers=false, label="", fillalpha=0.7, grid=false; color)
        if method == "SVD"
            @unpack m_difs, λs, rs = load(f_SVD_only_atm)
            iλ         = findfirst(λs .== λ0)
            ir         = findfirst(rs .== r0)
            boxplot!(p_box, m_difs[iλ, ir], outliers=false, ylabel=L"Cross-validation error $\epsilon\quad\mathrm{(m)}$", xticks=([1,2,3],["GP\n(only ATM)", "SVD\n(all)", "SVD\n(only ATM)"]), label="", fillalpha=0.7, grid=false; color)
        end
    end
    p_std = plot(p_std, legend=:topright, legend_foreground_color=nothing)
    GrIS1980s_DEM.panel_annotate!(p_std, l1)
    p_mean = plot(p_mean, legend=false, legend_foreground_color=nothing)
    p_box = plot(p_box, legend=:topright, legend_foreground_color=nothing)
    GrIS1980s_DEM.panel_annotate!(p_box, l3)
    return p_std, p_mean, p_box
end
p_std1, p_mean1, p_box1 = make_validation_plots("a", "b", "c", standardized=false)
p_std2, p_mean2, p_box2 = make_validation_plots("d", "e", "f", standardized=true)
p_both_lin = plot(p_std1, p_mean1, p_box1, p_std2, p_mean2, size=(2500,1400), left_margin=40mm, right_margin=15mm, bottom_margin=20mm, top_margin=20mm, dpi=300, layout=grid(2, 3, widths=(0.35, 0.35, 0.3)))
savefig(joinpath(fig_dir_main, "f04.png"))


##################################
# Uncertainty field GP, Figure 5 #
##################################

scalefontsizes()
scalefontsizes(1.3)
attr = (; grid=false, aspect_ratio=1, size=(500,900), xlims=(-7e5,8e5).*1e-3, ylims=(-3.32e6, -0.78e6).*1e-3, xaxis=false, yaxis=false, margin=10mm)
std_GP_std = NCDataset("output/reconstructions/interpolated_dh_std_GP.nc")["sigma_std"][:,:]
p1 = heatmap(x.*1e-3, y.*1e-3, std_GP_std', cmap=cgrad(:lapaz,rev=true), colorbar_title="\n(-)"; attr...)
GrIS1980s_DEM.panel_annotate!(p1, "a")
annotate!(p1, (xlims(p1)[1]+(xlims(p1)[2]-xlims(p1)[1])*0.55, ylims(p1)[1]+(ylims(p1)[2]-ylims(p1)[1])*1.07, text("STD of GP output \n(standardized)", 16)))
plot!(p1, outl, fillalpha=0, linewidth=0.4)
# add scalebar and North arrow
plot!(p1, [-7e2,-5e2], [-3e3, -3e3], color=:black, lw=3, label="")
annotate!(p1, [(-6e2, -2.85e3, text("200 km", :center, :Helvetica))])
plot!(p1, [-6e2,-6e2], [-2.6e3, -2.2e3], arrow=true, color=:black, lw=1, label="")
annotate!(p1, [(-6.4e2, -2.45e3, text("N", :right, :Helvetica))])
std_GP = NCDataset(get_rec_file_GP(grd))["std_uncertainty"][:,:]
p2 = heatmap(x.*1e-3, y.*1e-3, std_GP', clims=(0,20), cmap=cgrad(:lapaz,rev=true), colorbar_title="\n(m)"; attr...)
GrIS1980s_DEM.panel_annotate!(p2, "b")
annotate!(p2, (xlims(p2)[1]+(xlims(p2)[2]-xlims(p2)[1])*0.55, ylims(p2)[1]+(ylims(p2)[2]-ylims(p2)[1])*1.07, text("uncertainty map \n (de-standardized)", 16)))
plot!(p2, outl, fillalpha=0, linewidth=0.4)
plot(p1,p2, layout=(1,2), size=(1200, 700), left_margin=0mm, right_margin=[-20mm 10mm], bottom_margin=0mm, top_margin=20mm)
savefig(joinpath(fig_dir_main, "f05.png"))





# ###################
# # zoomed-in modes #
# ###################

# # Sermeq
# yb = -2310000
# xl = -227000
# # Helheim
# yb = -2600000
# xl =  260000
# ix = findall(xl.-10000 .< x .< xl+80000)
# iy = findall(yb.-8000 .< y .< yb+80000)
# slope_plot = zeros(size(dif))
# slope_plot[idif] .= bfield_1_m[idif]
# slope_plot[slope_plot .== 0] .= NaN
# heatmap(x[ix], y[iy], slope_plot[ix,iy]')

# d = load("output/reconstructions/SVD_components_g600_rec_without_W.jld2")
# # d = load(joinpath("output", "reconstructions", "SVD_components_g$(grd)_rec.jld2"))
# U = d["U"]

# b = zeros(length(x),length(y))
# ps = []
# for k in [1,5,10,20,50,80,100,200,300,600,700,750]
#     b[I_no_ocean] = U[:,k]
#     b[b .== 0] .= NaN
#     p_i = heatmap(x[ix], y[iy], b[ix,iy]', title="\n\n $k", cmap=:vik, clims=(-1e-2,1e-2), aspect_ratio=1, size=(900,800), margin=-20mm)
#     push!(ps, p_i)
# end
# plot(ps..., layout=(3,4), size=(2300,2000))
# savefig("modes_ii_Helheim.png")



##########################################
# SVD with vs without weights, Figure S6 #
##########################################


fs_SVD = [get_cv_file_SVD(grd, 70, only_atm=false),"output/validation/SVD_with_W_onlyatm_false.jld2"]
p_mean = plot(wsize=(GrIS1980s_DEM.wwidth, GrIS1980s_DEM.wheight))
p_std  = plot(wsize=(GrIS1980s_DEM.wwidth, GrIS1980s_DEM.wheight))
p_box  = plot()
cols   = palette(:batlow10)[[2,7]]
labels = ["no weights", "with weights"]
xlabel = "Elevation of reference DEM (m)"
# xticks = ([0.0, 50.,100.,150.], string.([0, 50,100,150]))
for (f_SVD, color, label) in zip(fs_SVD, cols, labels)
    @unpack idx, m_difs, λs, rs, nfiles, norms_UΣ, Σ = load(f_SVD)
    iλ = findfirst(λs .== λ0)
    ir = findfirst(rs .== 500)
    difs = m_difs[iλ, ir]
    dh_binned, bin_centers = GrIS1980s_DEM.bin_equal_bin_size(h_ref_m[idx], difs, 14)
    plot!(p_mean, bin_centers, mean.(dh_binned), ylabel=L"Error mean $\mu_\epsilon$ (m)", ls=:dot, marker=:circle; label, color, xlabel, attr...)
    hline!(p_mean, [0.0], color="grey", lw=3, z_order=1, label="", ls=:dash)
    plot!(p_std,    bin_centers, std.(dh_binned), ylabel=L"Error STD $\sigma_\epsilon\quad\mathrm{(m)}$", ls=:dot, marker=:circle; label, color, xlabel, attr...)
    boxplot!(p_box, difs, outliers=false, ylabel=L"Cross-validation error $\epsilon\quad\mathrm{(m)}$", xticks=([1,2],labels), label="", fillalpha=0.5, grid=false; color)
    if label == "with weights"
        global bin_centers = bin_centers
        global n_samples   = length.(dh_binned)
    end
end
p_std = plot(p_std, legend_foreground_color=nothing, legend=:topright, grid=false)
GrIS1980s_DEM.panel_annotate!(p_std, "a")
p_mean = plot(p_mean, legend=false, grid=false)
GrIS1980s_DEM.panel_annotate!(p_mean, "b")
p_tw = twinx(p_std)
bar!(p_tw, bin_centers, n_samples, color=:slategray, linecolor=:slategray, alpha=0.2, label="Relative bin size", grid=false, legend=:bottomleft, legend_foreground_color=:white, yaxis=false, right_margin=-30mm)
plot(p_std, p_mean, p_box ,size=(2000,700), left_margin=15mm, bottom_margin=15mm, top_margin=15mm, dpi=300, layout=grid(1, 3, widths=(0.38, 0.37, 0.25)))
savefig(joinpath(fig_dir_main, "fS06.png"))



##############################################################
# SVD error for different λ&r and singular values, Figure S5 #
##############################################################

# get file
f_SVD = get_cv_file_SVD(grd, 70, only_atm=false)
# plot absolute mean of cross-validation error for different λ and r values
attr = (;size=(900,700), margin=15Plots.mm, markersize=6, lw=3.5, markerstrokewidth=0)
p_mean = plot(xlabel="Elevation of reference DEM (m)", ylabel=L"Mean $\epsilon$ (m)")
p_std  = plot(xlabel="Elevation of reference DEM (m)", ylabel=L"Std $\epsilon$ (m)")
@unpack m_difs, λs, rs, Σ = load(f_SVD)
p = plot(xscale=:log10, xlabel="λ", xticks=10 .^(5:11), ylabel="Mean absolute error (m)", title="\n Cross-validation error", size=(1300,800), grid=false,
         legend=(0.15,0.92), palette = :tol_light, legend_foreground_color=nothing)
for i in eachindex(rs)
    md_abs = [abs.(md) for md in m_difs[:,i] ]
    plot!(p, λs, mean.(md_abs), label="r="*string(rs[i]), marker=:circle; attr...)
end
GrIS1980s_DEM.panel_annotate_xlog!(p, "a")
# plot singular values
p_sigm = plot(Σ.^2, yscale=:log10, yticks=[1e5, 1e7, 1e9, 1e11, 1e13, 1e15, 1e17], color=:black, ylabel=L"$\sigma_i^2$", xlabel=L"Mode index $i$", title="\n Singular values", label=""; attr...)
GrIS1980s_DEM.panel_annotate_ylog!(p_sigm, "b")
# save both in one figure
plot(p, p_sigm, size=(1800,600), dpi=300)
savefig(joinpath(fig_dir_main,"fS05.png"))


####################################################################
# SVD error for different number of training data files, Figure S7 #
####################################################################

# get file
fs_SVD = glob(get_cv_file_SVD(grd, "*", only_atm=false))
# plot
attr = (;margin=10Plots.mm, size=(GrIS1980s_DEM.wwidth,GrIS1980s_DEM.wheight), lw=1.8, markerstrokewidth=0, marker=:circle, ls=:dot, markersize=6, markeralpha=1.0)
p_mean = plot(wsize=(GrIS1980s_DEM.wwidth, GrIS1980s_DEM.wheight))
p_std  = plot(wsize=(GrIS1980s_DEM.wwidth, GrIS1980s_DEM.wheight))
cols   = palette(:batlow10)[1:2:end]
xlabel = "Elevation of reference DEM (m)"
for (f_SVD, color) in zip(fs_SVD, cols)
    @unpack idx, m_difs, λs, rs, nfiles, Σ = load(f_SVD)
    # m_B = zeros(size(h_ref_m))
    # m_B[I_no_ocean] .= norms_UΣ ./ sqrt(length(Σ)-1)
    if nfiles == 10   # doesn't have 500 modes
        iλ = findfirst(λs .== 1e6)
        ir = findfirst(rs .== 300)
    else
        iλ = findfirst(λs .== λ0)
        ir = findfirst(rs .== r0)
    end
    difs = m_difs[iλ, ir]
    # ii = findall(m_B[idx] .< 200)
    # dh_binned, bin_centers = GrIS1980s_DEM.bin_equal_bin_size(m_B[idx[ii]], difs[ii], 12)
    global dh_binned, bin_centers = GrIS1980s_DEM.bin_equal_bin_size(h_ref_m[idx], difs, 14)
    plot!(p_mean, bin_centers, mean.(dh_binned), legend_title=L"$m$", ylabel=L"Error mean $\mu_\epsilon$ (m)"; color, xlabel, attr...)
    hline!(p_mean, [0.0], color="grey", lw=3, z_order=1, label="", ls=:dash)
    plot!(p_std,    bin_centers, std.(dh_binned), legend_title=L"$m$", label=" $(length(Σ))", ylabel=L"Error STD $\sigma_\epsilon\quad\mathrm{(m)}$"; color, xlabel, attr...)
end
p_std = plot(p_std, legend_foreground_color=nothing, legend=:topright)
GrIS1980s_DEM.panel_annotate!(p_std, "a")
p_tw = twinx(p_std)
bar!(p_tw, bin_centers, length.(dh_binned), color=:slategray, linecolor=:slategray, alpha=0.2, label="Relative bin size", grid=false, legend=:bottomleft, legend_foreground_color=:white, yaxis=false, right_margin=-30mm)
p_mean = plot(p_mean, legend=false)
GrIS1980s_DEM.panel_annotate!(p_mean, "b")
plot(p_std, p_mean, size=(2000,700), left_margin=15mm, bottom_margin=15mm, top_margin=15mm, dpi=300)
savefig(joinpath(fig_dir_main, "fS07.png"))
