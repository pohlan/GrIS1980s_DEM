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
    f_dict = [get_cv_file_GP(grd, only_atm=standardized), get_cv_file_SVD(grd, 70, only_atm=standardized)]
    p_mean   = plot()
    p_std    = plot()
    p_box    = plot()
    for f in f_dict
        @unpack idx, grd, method = load(f)
        if method == "GP"
            @unpack difs, dists, binfield1, h_ref = load(f)
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
                difs       = standardize(difs, bfield_1_m[idx], h_ref_m[idx], subtract_mean=false)
                bin_metric = vcat(dists...)
            else
                bin_metric = h_ref_m[idx]
            end
        end
        # bin error vs elevation
        dh_binned, bin_centers = GrIS1980s_DEM.bin_equal_bin_size(bin_metric, difs, 12)
        # color
        color = standardized ? GrIS1980s_DEM.palette_dict[method*"_std"] : GrIS1980s_DEM.palette_dict[method]
        # plot mean and std errors
        label = standardized ? split(method,"_")[1]*" (only ATM)" : split(method,"_")[1]
        ylabel_mean = standardized ? L"Error mean $\mu_\epsilon$ (standardized)" : L"Error mean $\mu_\epsilon\quad\mathrm{(m)}$"
        xlabel = standardized ? "Distance to closest observation (km)" : "Elevation of reference DEM (m)"
        xscl    = standardized ? 1e-3 : 1
        plot!(p_mean, bin_centers.*xscl, mean.(dh_binned), ylabel=ylabel_mean, ls=:dot; color, xlabel, label, attr...) # y_foreground_color_text=color, y_guidefontcolor=color
        hline!(p_mean, [0.0], color="grey", lw=3, z_order=1, label="", ls=:dash) #, alpha=0.7; color)
        if !standardized
            plot!(p_std, bin_centers.*xscl, std.(dh_binned), ylabel=L"Error STD $\sigma_\epsilon\quad\mathrm{(m)}$", ls=:dot; color, xlabel, attr...)
        elseif method == "GP"
            plot!(p_std, bin_centers.*xscl, std.(dh_binned), ylabel=L"Error STD $\sigma_\epsilon$ (standardized)", ls=:dot, y_foreground_color_text=color, y_guidefontcolor=color; color, xlabel, attr...)
        else
            p_tw = twinx(p_std)
            plot!(p_tw, bin_centers.*xscl, std.(dh_binned), ylabel=L"Error STD $\sigma_\epsilon$ (standardized)", ls=:dot, y_foreground_color_text=color, y_guidefontcolor=color; color, attr...)
        end
        boxplot!(p_box, difs, outliers=false, ylabel=L"Cross-validation error $\epsilon\quad\mathrm{(m)}$", xticks=([1,2],["GP", "SVD"]), label="", fillalpha=0.7, grid=false; color)
        # uncertainty B-spline fit
        # insert!(bin_centers, 1, 0.0)
        # insert!(dh_binned, 1, [0.0,0.0])
        # sitp, rec_errors = uncertainty_from_cv(dh_binned, bin_centers, min_dists)
        # x_plot = 0.0:(diff(bin_centers)[1]*0.1):100e3
        # plot!(p_std, x_plot.*1e-3, sitp.(x_plot), z_order=1, label=split(method,"_")[1]*" B-spline fit"; color)
    end
    p_mean = plot(p_mean, legend_foreground_color=nothing)
    GrIS1980s_DEM.panel_annotate!(p_mean, l1)
    p_std = plot(p_std, legend=false, legend_foreground_color=nothing)
    GrIS1980s_DEM.panel_annotate!(p_std, l2)
    p_box = plot(p_box, legend=:topright, legend_foreground_color=nothing)
    GrIS1980s_DEM.panel_annotate!(p_box, l3)
    return p_mean, p_std, p_box
end
p_mean1, p_std1, p_box1 = make_validation_plots("a", "b", "c", standardized=false)
p_mean2, p_std2, p_box2 = make_validation_plots("d", "e", "f", standardized=true)
p_both_lin = plot(p_mean1, p_std1, p_box1, p_mean2, p_std2, p_box2, size=(2500,1400), left_margin=40mm, right_margin=15mm, bottom_margin=20mm, top_margin=20mm, dpi=300, layout=grid(2, 3, widths=(0.4, 0.4, 0.2)))
savefig(joinpath(fig_dir_main, "f04.png"))


##################################################
# Kriging errors for different maxns, Figure S1  #
##################################################



###############################
# SVD with vs without weights #
###############################


fs_SVD = [get_cv_file_SVD(grd, 70, only_atm=false),"output/validation/SVD_with_W_onlyatm_false.jld2"]
p_mean = plot(wsize=(GrIS1980s_DEM.wwidth, GrIS1980s_DEM.wheight))
p_std  = plot(wsize=(GrIS1980s_DEM.wwidth, GrIS1980s_DEM.wheight))
p_box  = plot()
cols   = palette(:batlow10)[[2,7]]
labels = ["no weights", "with weights"]
xlabel = "Elevation of reference DEM (m)"
xticks = ([0.0, 50.,100.,150.], string.([0, 50,100,150]))
for (f_SVD, color, label) in zip(fs_SVD, cols, labels)
    @unpack idx, m_difs, λs, rs, nfiles, norms_UΣ, Σ = load(f_SVD)
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
    # ii = findall(m_B[idx] .< 150)
    # dh_binned, bin_centers = GrIS1980s_DEM.bin_equal_bin_size(m_B[idx[ii]], difs[ii], 12)
    dh_binned, bin_centers = GrIS1980s_DEM.bin_equal_bin_size(h_ref_m[idx], difs, 14)
    plot!(p_mean, bin_centers, mean.(dh_binned), ylabel=L"Mean $\epsilon$ (m)", ls=:dot; label, color, xlabel, xticks, attr...)
    hline!(p_mean, [0.0], color="grey", lw=3, z_order=1, label="", ls=:dash)
    plot!(p_std,    bin_centers, std.(dh_binned), ylabel=L"Error standard deviation $\sigma_\epsilon\quad\mathrm{(m)}$", ls=:dot; label, color, xlabel, xticks, attr...)
    boxplot!(p_box, difs, outliers=false, ylabel=L"Cross-validation error $\epsilon\quad\mathrm{(m)}$", xticks=([1,2],labels), label="", fillalpha=0.5, grid=false; color)
    if label == "with weights"
        global bin_centers = bin_centers
        global n_samples   = length.(dh_binned)
    end
end
p_mean = plot(p_mean, legend_foreground_color=nothing, legend=:topright, grid=false)
GrIS1980s_DEM.panel_annotate!(p_mean, "a")
p_std = plot(p_std, legend=false, grid=false)
GrIS1980s_DEM.panel_annotate!(p_std, "b")
# p_tw = twinx(p_std)
# bar!(p_tw, bin_centers, n_samples, color=:slategray, linecolor=:slategray, alpha=0.2, label="relative sample size", grid=false, legend=:top, legend_foreground_color=nothing, yaxis=false, right_margin=-30mm)
plot(p_mean, p_std, p_box ,size=(2000,700), margin=15Plots.mm, dpi=300, layout=grid(1, 3, widths=(0.38, 0.37, 0.25)))
savefig(joinpath(fig_dir_main, "fS0xx.png"))



##############################################################
# SVD error for different λ&r and singular values, Figure S3 #
##############################################################

# get file
f_SVD = get_cv_file_SVD(grd, 70, only_atm=false)
# f_SVD = "output/validation/cv_file_it3.jld2"
# f_SVD = "output/validation/cv_1e5.3_gr600_SVD_nfiles10.jld2"
# f_SVD = "output/validation/SVD_with_W_onlyatm_false.jld2"
# plot absolute mean of cross-validation error for different λ and r values
attr = (;size=(900,700), margin=15Plots.mm, markersize=6, lw=3.5, markerstrokewidth=0)
p_mean = plot(xlabel="Elevation of reference DEM (m)", ylabel=L"Mean $\epsilon$ (m)")
p_std  = plot(xlabel="Elevation of reference DEM (m)", ylabel=L"Std $\epsilon$ (m)")
@unpack m_difs, λs, rs, Σ = load(f_SVD)
p = plot(xscale=:log10, xlabel="λ", xticks=10 .^(5:11), ylabel="Mean absolute error (m)", title="\n Cross-validation error", size=(1300,800), grid=false,
         legend=(0.15,0.75), palette = :tol_light, legend_foreground_color=nothing)
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
savefig(joinpath(fig_dir_main,"fS03.png"))


####################################################################
# SVD error for different number of training data files, Figure S4 #
####################################################################

# get file
fs_SVD = glob(get_cv_file_SVD(grd, "*", only_atm=false))
# plot
attr = (;margin=10Plots.mm, size=(GrIS1980s_DEM.wwidth,GrIS1980s_DEM.wheight), lw=1.8, markerstrokewidth=0, marker=:circle, ls=:dot, markersize=6, markeralpha=1.0)
p_mean = plot(wsize=(GrIS1980s_DEM.wwidth, GrIS1980s_DEM.wheight))
p_std  = plot(wsize=(GrIS1980s_DEM.wwidth, GrIS1980s_DEM.wheight))
cols   = palette(:batlow10)[1:2:end]
xlabel = "Standard deviation of model realizations (m)"
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
    dh_binned, bin_centers = GrIS1980s_DEM.bin_equal_bin_size(h_ref_m[idx], difs, 14)
    plot!(p_mean, bin_centers, mean.(dh_binned), legend_title=L"$m$", label=" $(length(Σ))", ylabel=L"Mean $\epsilon$ (m)"; color, xlabel, attr...)
    hline!(p_mean, [0.0], color="grey", lw=3, z_order=1, label="", ls=:dash)
    plot!(p_std,    bin_centers, std.(dh_binned), legend_title=L"$m$", ylabel=L"Error standard deviation $\sigma_\epsilon\quad\mathrm{(m)}$"; color, xlabel, attr...)
end
p_mean = plot(p_mean, legend_foreground_color=nothing, legend=:bottomright)
GrIS1980s_DEM.panel_annotate!(p_mean, "a")
p_std = plot(p_std, legend=false)
GrIS1980s_DEM.panel_annotate!(p_std, "b")
plot(p_mean, p_std,size=(2000,700), margin=15Plots.mm, dpi=300)
savefig(joinpath(fig_dir_main, "fS04.png"))
