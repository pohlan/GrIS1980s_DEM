# define output path for other figures not in paper
fig_dir_others = joinpath("output", "data_preprocessing", "figures")
mkpath(fig_dir_others)

##############################
# Standardization, Figure 1  #
##############################

# plotting parameters
Plots.scalefontsizes()
Plots.scalefontsizes(1.7)
cmap_div  = :RdBu
cmap_else = :batlow
wheight   = 700
wwidth    = 1000

# std
heatmap(bin_centers_1, bin_centers_2, nmads', cmap=:batlow)
savefig(joinpath(fig_dir_others, "nmads_2Dbinning.png"))

# bias
heatmap(bin_centers_1, bin_centers_2, meds')
savefig(joinpath(fig_dir_others, "medians_2Dbinning.png"))

# qqplot
plot(
    StatsPlots.qqplot(StatsPlots.Normal(), df_all.dh_detrend, title="standardized with binning", ylims=(-8,8)),
    StatsPlots.qqplot(StatsPlots.Normal(), (df_all.dh .- mean(df_all.dh))./std(df_all.dh), title="standardized without binning", ylims=(-8,8)), size=(1200,700)
    )
savefig(joinpath(fig_dir_others,"qqplot.png"))

# interpolations
x1 = range(bin_centers_1[1], bin_centers_1[end], length=10000)
x2 = range(bin_centers_2[1], bin_centers_2[end], length=1000)
itp_var = GrIS1980s_DEM.get_itp_interp(bin_centers_1, bin_centers_2, nmads)
p_std  = heatmap(x1, x2, itp_var.(x1, x2')',  xlabel="Slope (°)", ylabel="Elevation (m)", colorbar_title=L"$\sigma_{\Delta h}\,\mathrm{(m)}$", wsize=(wwidth, wheight), margin=12Plots.mm, cmap=cmap_else)
GrIS1980s_DEM.panel_annotate!(p_std, "a")
itp_bias = GrIS1980s_DEM.get_itp_interp(bin_centers_1, bin_centers_2, meds)
p_bias = heatmap(x1, x2, itp_bias.(x1, x2')', xlabel="Slope (°)", ylabel="Elevation (m)", colorbar_title=L"$\overline{\Delta h}\,\mathrm{(m)}$", wsize=(wwidth, wheight), margin=12Plots.mm, cmap=cmap_div, clims=(-40,40))
GrIS1980s_DEM.panel_annotate!(p_bias, "b")

# histogram
p_hist = histogram(df_all.dh_detrend, label="Standardized \nobservations", xlims=(-7,7), yticks=false, xlabel="Elevation difference (m)", ylabel="Frequency", normalize=:pdf, nbins=800, color=:cornflowerblue, wsize=(wwidth, wheight), linecolor=:cornflowerblue, margin=12Plots.mm)
plot!(p_hist, Normal(), lw=2.5, label="Normal distribution", color="black", foreground_color_legend = nothing)
GrIS1980s_DEM.panel_annotate!(p_hist, "c")

# variogram
xvals, yvals = values(gamma)
varg = GrIS1980s_DEM.get_var(gamma; adjust_sill=false)
p_varg = scatter(ustrip.(xvals) .* 1e-3, yvals ./ sill(varg), label="Empirical variogram", color=:cornflowerblue, markerstrokewidth=0, wsize=(wwidth,wheight), xlabel="Distance (km)", ylabel=L"\gamma\,\,\mathrm{(m^2)}", margin=12Plots.mm)
plot!(p_varg, [1e-5; ustrip.(xvals)] .* 1e-3, varg.([1e-5; ustrip.(xvals)]) ./ sill(varg), label="Variogram fit", lw=3.5, ylims=(0,1.3), color=:black, foreground_color_legend=nothing)
GrIS1980s_DEM.panel_annotate!(p_varg, "d")

# FIGURE 1
p_all = plot(p_std, p_bias, p_hist, p_varg, wsize=(2500, 1200), left_margin=30Plots.mm, dpi=300)
savefig(joinpath(fig_dir_main, "f01.png"))


##############################
# Heatmaps of data, Figure 2 #
##############################

function heatmap_from_df(df_all, sm::Symbol, x, y, dims::Tuple, fname; clims=(-4,4), title)
    # plotting parameters
    Plots.scalefontsizes()
    Plots.scalefontsizes(1.3)
    wsize          = (1000,1200)
    xlabel         = "Easting (m)"
    ylabel         = "Northing (m)"
    colorbar_title = "(m)"
    cmap           = :RdBu
    # plot raster with aerodem data
    id_df_aero = findall(df_all.source .== :aerodem .|| df_all.source .== "aerodem")
    m_plot = zeros(dims)
    m_plot[df_all.idx[id_df_aero]] .= df_all[!,sm][id_df_aero]
    i_nans = m_plot .== 0
    i_nans[I_no_ocean] .= false
    m_plot[i_nans] .= NaN
    heatmap(x, y, m_plot'; cmap, dpi=300)
    # scatter atm point data
    id_df_atm = findall(df_all.source .== :atm .|| df_all.source .== "atm")
    scatter!(df_all.x[id_df_atm], df_all.y[id_df_atm], marker_z=df_all[!,sm][id_df_atm], label="", markersize=1.2, markerstrokewidth=0, aspect_ratio=1, xlims=(-7e5,8e5), ylims=(-3.32e6, -0.78e6), grid=false, bottommargin=0Plots.cm, rightmargin=(8/5)Plots.cm, leftmargin=(10/5)Plots.cm; wsize, xlabel, ylabel, colorbar_title, title, clims, cmap)
    plot!(outl, fill=nothing, lw=0.2)
    savefig(fname)
    return
end

# plot before standardizing
heatmap_from_df(df_all, :dh, x, y, size(h_ref_m), joinpath(fig_dir_others,"data_non-standardized.png"), clims=(-20,20), title=L"Non-standardized observations $\Delta h$")

# plot after standardizing (FIGURE 2)
heatmap_from_df(df_all, :dh_detrend, x, y, size(h_ref_m), joinpath(fig_dir_main,"f02.png"), clims=(-3.0,3.0), title=L"Standardized observations $\Delta h_\mathrm{std}$")
