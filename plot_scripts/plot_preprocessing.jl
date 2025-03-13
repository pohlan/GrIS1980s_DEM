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
Plots.heatmap(bin_centers_1, bin_centers_2, nmads', cmap=:batlow)
Plots.savefig(joinpath(fig_dir_others, "nmads_2Dbinning.png"))

# bias
Plots.heatmap(bin_centers_1, bin_centers_2, meds')
Plots.savefig(joinpath(fig_dir_others, "medians_2Dbinning.png"))

# qqplot
Plots.plot(
    StatsPlots.qqplot(StatsPlots.Normal(), df_all.dh_detrend, title="standardized with binning", ylims=(-8,8)),
    StatsPlots.qqplot(StatsPlots.Normal(), (df_all.dh .- mean(df_all.dh))./std(df_all.dh), title="standardized without binning", ylims=(-8,8)), size=(1200,700)
    )
Plots.savefig(joinpath(fig_dir_others,"qqplot.png"))

# interpolations
x1 = range(bin_centers_1[1], bin_centers_1[end], length=10000)
x2 = range(bin_centers_2[1], bin_centers_2[end], length=1000)
itp_var = svd_IceSheetDEM.get_itp_interp(bin_centers_1, bin_centers_2, nmads)
p_std  = heatmap(x1, x2, itp_var.(x1, x2')',  xlabel="Slope (°)", ylabel="Elevation (m)", colorbar_title=L"$\sigma_{\Delta h}\,\mathrm{(m)}$", wsize=(wwidth, wheight), margin=12Plots.mm, cmap=cmap_else)
annotate!(p_std, (xlims(p_std)[1], ylims(p_std)[2]*1.05, text(L"\textbf{a}", :left, 27)))  #  "Times Bold", "Helvetica Bold"
itp_bias = svd_IceSheetDEM.get_itp_interp(bin_centers_1, bin_centers_2, meds)
p_bias = heatmap(x1, x2, itp_bias.(x1, x2')', xlabel="Slope (°)", ylabel="Elevation (m)", colorbar_title=L"$\overline{\Delta h}\,\mathrm{(m)}$", wsize=(wwidth, wheight), margin=12Plots.mm, cmap=cmap_div, clims=(-40,40))
annotate!(p_bias, (xlims(p_bias)[1], ylims(p_bias)[2]*1.05, text(L"\textbf{b}", :left, 27)))  #  "Times Bold", "Helvetica Bold"

# histogram
p_hist = histogram(df_all.dh_detrend, label="Standardized \nobservations", xlims=(-7,7), yticks=false, xlabel="Elevation difference (m)", ylabel="Frequency", normalize=:pdf, nbins=1500, color=:cornflowerblue, wsize=(wwidth, wheight), linecolor=nothing, margin=12Plots.mm)
plot!(p_hist, Normal(), lw=2.5, label="Normal distribution", color="black", foreground_color_legend = nothing)
annotate!(p_hist, (xlims(p_hist)[1], ylims(p_hist)[2]*1.05, text(L"\textbf{c}", :left, 27)))  #  "Times Bold", "Helvetica Bold"

# variogram
xvals, yvals = values(gamma)
varg = svd_IceSheetDEM.get_var(gamma; adjust_sill=false)
p_varg = scatter(ustrip.(xvals) .* 1e-3, yvals ./ sill(varg), label="Empirical variogram", color=:cornflowerblue, markerstrokewidth=0, wsize=(wwidth,wheight), xlabel="Distance (km)", ylabel=L"\gamma\,\,\mathrm{(m^2)}", margin=12Plots.mm)
plot!(p_varg, [1e-5; ustrip.(xvals)] .* 1e-3, varg.([1e-5; ustrip.(xvals)]) ./ sill(varg), label="Variogram fit", lw=3.5, ylims=(0,1.3), color=:black, foreground_color_legend=nothing)
annotate!(p_varg, (xlims(p_varg)[1], ylims(p_varg)[2]*1.05, text(L"\textbf{d}", :left, 27)))  #  "Times Bold", "Helvetica Bold"

# FIGURE 1
p_all = plot(p_std, p_bias, p_hist, p_varg, wsize=(2500, 1200), left_margin=30Plots.mm)
Plots.savefig(joinpath(fig_dir_main, "Figure1.png"))


##############################
# Heatmaps of data, Figure 2 #
##############################

function heatmap_from_df(df_all, sm::Symbol, x, y, dims::Tuple, fname; clims=(-4,4), title)
    # plotting parameters
    Plots.scalefontsizes()
    Plots.scalefontsizes(6.5)
    wsize          = (5000,6000)
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
    heatmap(x, y, m_plot'; cmap) #, aspect_ratio=1, xlims=(-7e5,8e5); clims)
    # scatter atm point data
    id_df_atm = findall(df_all.source .== :atm .|| df_all.source .== "atm")
    scatter!(df_all.x[id_df_atm], df_all.y[id_df_atm], marker_z=df_all[!,sm][id_df_atm], label="", markersize=6.0, markerstrokewidth=0, aspect_ratio=1, xlims=(-7e5,8e5), ylims=(-3.32e6, -0.78e6), grid=false, bottommargin=0Plots.cm, rightmargin=8Plots.cm, leftmargin=10Plots.cm; wsize, xlabel, ylabel, colorbar_title, title, clims, cmap)
    plot!(shp, fill=nothing, lw=1)
    savefig(fname)
    return
end

# plot before standardizing
heatmap_from_df(df_all, :dh, x, y, size(h_ref), joinpath(fig_dir_others,"data_non-standardized.png"), clims=(-20,20), title=L"Non-standardized observations $\Delta h$")

# plot after standardizing (FIGURE 2)
heatmap_from_df(df_all, :dh_detrend, x, y, size(h_ref), joinpath(fig_dir_main,"Figure2.png"), clims=(-3.0,3.0), title=L"Standardized observations $\Delta h_\mathrm{std}$")
