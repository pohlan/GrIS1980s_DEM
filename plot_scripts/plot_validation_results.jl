using JLD2, UnPack, Glob, Plots, LaTeXStrings

# plot
Plots.scalefontsizes()
Plots.scalefontsizes(1.9)
attr = (;margin=10Plots.mm, size=(900,700), lw=2, markerstrokewidth=0.5, marker=:circle)

output_dir = "output/validation/"
fig_dir    = joinpath(output_dir, "figures/")
mkpath(fig_dir)

f_dict = glob("output/validation/dict_cv_block_1e5.3_gr1200_*.jld2")

# choose SVD parameterization (λ and r)
i_SVD = 3

p_median = plot()
p_std    = plot()
p_variog = plot()
for f in f_dict
    @unpack xc, yc, idx, h_ref, gr, method = load(f)
    if method == "SVD"
        @unpack m_difs, λs, rs, dict = load(f)
        ## for SVD: plot mean, norm etc for different λ and r values
        # see https://juliagraphics.github.io/ColorSchemes.jl/stable/basics/#Pre-defined-schemes-1 for more pre-defined, color-blind friendly schemes
        for m in keys(dict)
            p = plot(xscale=:log10, xlabel="λ", ylabel=m, size=(1300,800), leftmargin=10Plots.mm, topmargin=10Plots.mm, bottommargin=10Plots.mm, legend=:top, palette = :tol_light)
            for i in eachindex(rs)
                plot!(λs, abs.(dict[m][:,i]), label="r="*string(rs[i]), marker=:circle, markersize=6, markerstrokewidth=0, lw=3.5)
            end
            plot(p)
            savefig(joinpath(fig_dir,m*"_abs_gr$(gr)_SVD.png"))
        end
        difs = m_difs[i_SVD]
    else
        @unpack difs = load(f)
    end

    # calculate de-standardized difs
    dif_destd = destandardize(difs, h_ref, add_mean=false)

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
    plot!(p_median, bin_centers, median.(dh_binned), label=method, xlabel="Elevation of reference DEM (m)", ylabel=L"Median $\epsilon$ (m)"; attr...)
    plot!(p_std,    bin_centers, std.(dh_binned), label=method, xlabel="Elevation of reference DEM (m)", ylabel="Std (m)"; attr...)

    # variogram
    data = svd_IceSheetDEM.make_geotable(dif_destd, xc, yc)
    U = data |> UniqueCoords()
    gamma = EmpiricalVariogram(U, :Z; estimator=:cressie, nlags=1000,  maxlag=1e6)
    id_plot = findall(gamma.ordinate .!= 0.0)
    scatter!(p_variog, gamma.abscissa[id_plot], gamma.ordinate[id_plot], markersize=2, label=method; attr...)
end
plot(p_median)
savefig(joinpath(fig_dir, "median_vs_elevation.png"))
plot(p_std)
savefig(joinpath(fig_dir, "std_vs_elevation.png"))
plot(p_variog)
savefig(joinpath(fig_dir, "error_variogram.png"))


### scatter plot to check correlation with h_kriging - h_SVD

# idx = vcat(idxs...)
# f_kriging = "output/geostats_interpolation/kriging/rec_kriging_g1200.nc"
# h_kriging = NCDataset(f_kriging)["surface"][:,:]
# f_SVD     = "output/SVD_reconstruction/rec_lambda_1e7_g1200_r200.nc"
# h_SVD = NCDataset(f_SVD)["surface"][:,:]
# dif_recs = h_kriging .- h_SVD
# xgrid = NCDataset(f_SVD)["x"][:,:]
# ygrid = NCDataset(f_SVD)["y"][:,:]

# itp = interpolate((xgrid, ygrid), dif_recs, Gridded(Linear()))
# df_atm = df_all[df_all.source .== "atm",:]
# dif_recs_pts = itp.(df_atm.x, df_atm.y)

# id_aero = findall(.!ismissing.(df_all.idx[idx]))
# drec = vcat(dif_recs[df_all.idx[idx[id_aero]]], dif_recs_pts)

# dif_cv = destandardize(difs, df_all.h_ref[idx], add_mean=false)
# histogram2d(drec, dif_cv, markerstrokewidth=0, markersize=0.5, cmap=:curl, label="", clims=(0,200), colorbar_ticks=false, colorbar_title="Density"; attr...)
# b = .!ismissing.(drec)
# yy = -280:280
# plot!(yy, yy, color="black", xlims=(-300,300), ylims=(-300,300), lw=3, xlabel=L"\delta h_\mathrm{kriging-SVD} (m)", ylabel=L"\epsilon_\mathrm{kriging} (m)", margin=10Plots.mm, size=(900,700), label="")
# cr = round(cor(drec[b], difs[b]), digits=2)
# title!("Correlation = $cr")
# savefig(joinpath(fig_dir, "correlation_error_dh_kriging.png"))
