# define output path for other figures not in paper
fig_dir_others = joinpath("output", "data_preprocessing", "figures")
mkpath(fig_dir_others)

##############################
# Standardization, Figure 2  #
##############################

# plotting parameters
Plots.scalefontsizes()
Plots.scalefontsizes(2.1)
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

# plot bin statistics: std, median and nsamples
# make polygon outlining the non-interpolated bins
bin_edges_1 = [bin_centers_1[1]-0.5*diff(bin_centers_1)[1]; bin_centers_1[2:end].-0.5*diff(bin_centers_1)]
bin_edges_2 = [bin_centers_2[1]-0.5*diff(bin_centers_2)[1]; bin_centers_2[2:end].-0.5*diff(bin_centers_2)]
xpts = bin_edges_1[[1,1,6,6,7,7,8,8,7,7,5,5,4,4,3,3,2,2,1,1,2,2,1,1,2,2,1]]
ypts = bin_edges_2[[3,1,1,2,2,8,8,15,15,19,19,21,21,20,20,24,24,26,26,14,14,11,11,10,10,3,3]]
pt_crds = [(xp,yp) for (xp,yp) in zip(xpts,ypts)]
df_pts   = DataFrame(geometry=AG.createpolygon(pt_crds))
pol = df_pts.geometry
# add nans to heatmap field
nmads_withnans = copy(nmads); nmads_withnans[nsamples_bins.==1] .= NaN
meds_withnans = copy(meds); meds_withnans[nsamples_bins.==1] .= NaN
ic = sum([all(isnan.(c)) for c in eachcol(nmads_withnans)])
ir = sum([all(isnan.(c)) for c in eachrow(nmads_withnans)])
function bin_plot(field_withnans, colorbar_title, cmap, panel_letter, yticks, lm; clims=nothing, padding=6)
    p_bins = plot(xlabel = L"Slope $s_\mathrm{ref}$ (°)", ylabel=L"Elevation $h_\mathrm{ref}$ (m)"; colorbar_title, right_margin=0Plots.mm, size=(wwidth,wheight))
    heatmap!(p_bins, bin_centers_1[1:end-ir], bin_centers_2[1:end-ic], field_withnans[1:end-ir,1:end-ic]'; cmap, colorbar=false, clims)
    plot!(p_bins, pol, lw=2, fillalpha=0, linecolor=:black)
    GrIS1980s_DEM.panel_annotate!(p_bins, panel_letter)
    # colorbar
    vmax = maximum(field_withnans[.!isnan.(field_withnans)])
    xx = isnothing(clims) ? range(0,vmax,500) : range(clims...,500)
    zz = zero(xx)' .+ xx
    p_c = heatmap(xx, xx, zz, ticks=false, ratio=25, legend=false, fc=cgrad(cmap), lims=extrema(xx), framestyle=:box, left_margin=-10Plots.mm, right_margin=50Plots.mm, top_margin=30Plots.mm, bottom_margin=20Plots.mm)
    ticktxt = string.(yticks)
    [annotate!(p_c, vmax*1.1, yi, text("- $ti", 20, "Helvetica", :left)) for (yi,ti) in zip(yticks,ticktxt)]
    annotate!(p_c, vmax*padding, mean(xx), text(colorbar_title, 23, rotation=90))
    # plot both together
    pboth = plot(p_bins, p_c, left_margin=10Plots.mm, right_margin=[-lm*mm -10mm])
return pboth
end
# plot for std, median and nsamples
p_bins_std = bin_plot(nmads_withnans, L"$\sigma_{\Delta h}\,\mathrm{(m)}$", cmap_else, "a", 10:10:60, 60, padding=6)
p_bins_med = bin_plot(meds_withnans, L"$\overline{\Delta h}\,\mathrm{(m)}$", cmap_div, "b", -30:10:30, clims=(-40,40), 60, padding=30)
nsamples_bins = Float64.(nsamples_bins)
nsamples_bins[nsamples_bins.==1] .= NaN
p_bins_density = bin_plot(nsamples_bins, "# samples per bin", :dense, "c", 2000:2000:10000, 60, padding=8)

# plot top panels together
p_top = plot(p_bins_std, p_bins_med, p_bins_density, size=(wwidth*3, wheight), layout=grid(1,3, widths=(3/9,3/9,3/9)), left_margin=[18mm 18mm 18mm]) #, right_margin=[-2mm -2mm -10mm]) #, size=(wwidth*1.5, wheight), left_margin=18Plots.mm, right_margin=18Plots.mm, top_margin=14Plots.mm, bottom_margin=12Plots.mm)

# histogram
p_hist = histogram(df_all.dh_detrend, label="Standardized \nobservations", xlims=(-7,7), yticks=false, xlabel="Elevation difference (-)", ylabel="Frequency", normalize=:pdf, nbins=800, color=:cornflowerblue, linecolor=:cornflowerblue)
plot!(p_hist, Normal(), lw=2.5, label="Normal distribution", color="black", foreground_color_legend = nothing)
GrIS1980s_DEM.panel_annotate!(p_hist, "d")

# variogram
xvals, yvals = values(gamma)
varg = GrIS1980s_DEM.get_var(gamma; adjust_sill=false)
p_varg1 = scatter(ustrip.(xvals) .* 1e-3, yvals ./ sill(varg), label="Empirical variogram", color=:cornflowerblue, markerstrokewidth=0, markersize=5, xlabel="Distance (km)", ylabel=L"\gamma", xlims=(0, 150))
plot!(p_varg1, [1e-5; ustrip.(xvals)] .* 1e-3, varg.([1e-5; ustrip.(xvals)]) ./ sill(varg), label="Variogram fit", lw=3.5, ylims=(0,1.3), color=:black, foreground_color_legend=nothing)
GrIS1980s_DEM.panel_annotate!(p_varg1, "e")
p_varg = p_varg1

p_bot = plot(p_hist, p_varg, right_margin=[40mm 10mm], left_margin=30mm)

# FIGURE 2
p_all = plot(p_top, p_bot, layout=(2,1), dpi=300, top_margin=20mm, bottom_margin=20mm, size=(wwidth*3, wheight*2))
savefig(joinpath(fig_dir_main, "f02.png"))

#############################################
# Variogram over stable terrain, Figure S4  #
#############################################
xvals, yvals = values(gamma_error)
varg = GrIS1980s_DEM.get_var(gamma_error, nVmax=1, adjust_sill=false)
p_varg = scatter(ustrip.(xvals) .* 1e-3, yvals ./ sill(varg), label="Empirical variogram", color=:darkorchid, markerstrokewidth=0, markersize=7, xlabel="Distance (km)", ylabel=L"\gamma", size=(GrIS1980s_DEM.wwidth,GrIS1980s_DEM.wheight), margin=6mm)
plot!(p_varg, [1e-5; ustrip.(xvals)] .* 1e-3, varg.([1e-5; ustrip.(xvals)]) ./ sill(varg), label="Variogram fit", lw=3.5, ylims=(0,1.3), color=:black, foreground_color_legend=nothing, legend_background_color=nothing, legend=:bottomright)
savefig(joinpath(fig_dir_main, "fS04.png"))


#############################################
# Interpolated/Extrapolated bins, Figure S2 #
#############################################

# fine-scale x values
x1 = range(bin_edges_1[1], bin_edges_1[end], length=10000)
x2 = range(bin_edges_2[1], bin_edges_2[end], length=1000)
# std
itp_var = GrIS1980s_DEM.get_itp_interp(bin_centers_1, bin_centers_2, nmads)
p_std  = heatmap(x1, x2, itp_var.(x1, x2')',  xlabel=L"Slope $s_\mathrm{ref}$ (°)", ylabel=L"Elevation $h_\mathrm{ref}$ (m)", colorbar_titlefontsize=25, colorbar_title=L"$\sigma_{\Delta h}\,\mathrm{(m)}$", wsize=(wwidth, wheight), cmap=cmap_else)
plot!(p_std, pol, lw=2, fillalpha=0, linecolor=:black)
GrIS1980s_DEM.panel_annotate!(p_std, "a")
# median
itp_bias = GrIS1980s_DEM.get_itp_interp(bin_centers_1, bin_centers_2, meds)
p_bias = heatmap(x1, x2, itp_bias.(x1, x2')', xlabel=L"Slope $s_\mathrm{ref}$ (°)", ylabel=L"Elevation $h_\mathrm{ref}$ (m)", colorbar_titlefontsize=25, colorbar_title=L"$\overline{\Delta h}\,\mathrm{(m)}$", wsize=(wwidth, wheight), cmap=cmap_div, clims=(-40,40))
plot!(p_bias, pol, lw=2, fillalpha=0, linecolor=:black)
GrIS1980s_DEM.panel_annotate!(p_bias, "b")
# subplots
plot(p_std, p_bias, size=(wwidth*2, wheight), margin=18mm)
savefig(joinpath(fig_dir_main, "fS02.png"))


####################################
# Variogram supplements, Figure S3 #
####################################


# aerodem vs ATM, separate variograms
gms_aero_atm = []
for (source,lab) in zip(["aerodem", "atm"],["AeroDEM", "ATM"])
    i_g = findall(df_all.source .== source)
    gm_i = GrIS1980s_DEM.emp_variogram(df_all.x[i_g], df_all.y[i_g], df_all.dh_detrend[i_g]; maxlag=2e5, nlags=200)
    push!(gms_aero_atm, gm_i)
end
# plot the two
palette = Plots.palette(Plots.palette(:Archambault)[[1,4]])
p_varg_atm_aero = plot(size=(wwidth,wheight), xlims=(-10,190), legend=:bottomright; palette)
for (gm_i, source) in zip(gms_aero_atm, ["AeroDEM", "ATM"])
    xv, yv = values(gm_i)
    ysc = yv[findmin(abs.(ustrip.(xv) .- 150e3))[2]]
    scatter!(p_varg_atm_aero, ustrip.(xv) .* 1e-3, yv ./ ysc, label=source, markersize=4, markerstrokewidth=0, xlabel="Distance (km)", ylabel=L"\gamma")
end
# add the ice-sheet wide variogram used for kriging (ATM + part of AeroDEM)
xv, yv = values(gamma)
varg = GrIS1980s_DEM.get_var(gamma, adjust_sill=false)
ysc = yv[findmin(abs.(ustrip.(xv) .- 150e3))[2]]
scatter!(p_varg_atm_aero, ustrip.(xv) .* 1e-3, yv ./ ysc, markersize=4, markerstrokewidth=0, label="AeroDEM + ATM", color=:black, foreground_color_legend=nothing)
plot!(p_varg_atm_aero, [1e-5; ustrip.(xv)] .* 1e-3, varg.([1e-5; ustrip.(xv)]) ./ varg(150e3), color=:black, label="Fitted model", lw=3.5)
GrIS1980s_DEM.panel_annotate!(p_varg_atm_aero, "a")
# map
aero_atm_map = plot(aspect_ratio=1, size=(wheight,wwidth), xlims=(-7e5,8e5).*1e-3, ylims=(-3.32e6, -0.78e6).*1e-3, xaxis=false, yaxis=false; palette)
for source in ["aerodem", "atm"]
    i_g = findall(df_all.source .== source)
    scatter!(aero_atm_map, df_all.x[i_g].*1e-3, df_all.y[i_g].*1e-3, markerstrokewidth=0, markersize=0.5, aspect_ratio=1, label="")
end
GrIS1980s_DEM.panel_annotate!(aero_atm_map, "b")
# scalebar and north arrow
plot!(aero_atm_map, [-7e2,-5e2], [-3e3, -3e3], color=:black, lw=3, label="")
annotate!(aero_atm_map, [(-6e2, -2.85e3, text("200 km", :center, :Helvetica, 20))])
plot!(aero_atm_map, [-6e2,-6e2], [-2.6e3, -2.2e3], arrow=true, color=:black, lw=1, label="")
annotate!(aero_atm_map, [(-6.4e2, -2.45e3, text("N", :right, :Helvetica, 20))])
p_atm_aero = plot(p_varg_atm_aero, aero_atm_map, layout=grid(1,2,widths=(0.65, 0.35)), size=(wwidth*2, wheight), right_margin=[0mm 10mm], left_margin=[20mm 5mm], bottom_margin=15mm)

# divide Greenland into sections
# take only a subsample of aerodem
i_aero = findall(df_all.source .== "aerodem")
i_atm  = findall(df_all.source .== "atm")
# n_cells_interior = length(I_no_ocean) - length(i_g)
# dens_atm = length(i_g2)/n_cells_interior # around 37 times higher density for aerodem
n_aero = round(Int, 0.2 * length(i_aero))
i_sampl = [rand(1:length(i_aero), n_aero); i_atm]
df_not_all = df_all[i_sampl,:]

ymid = -2e6
is = [findall(df_not_all.y .<= ymid),
      findall(ymid .< df_not_all.y)]
xvals_4 = []
yvals_4 = []
# rs_aero = [0.5,0.5,0.2]
for i_g in is
    gm = GrIS1980s_DEM.emp_variogram(df_not_all.x[i_g], df_not_all.y[i_g], df_not_all.dh_detrend[i_g]; maxlag=1.8e5, nlags=200)
    xvals, yvals = values(gm)
    push!(xvals_4, xvals); push!(yvals_4, yvals)
end
# plot the sections
palette = Plots.palette(Plots.palette(:tableau_temperature)[[3,6]])
p4_varg = plot(legend=:topleft, size=(wwidth,wheight), xlims=(-3,150), ylims=(0,1.2), foreground_color_legend=nothing; palette) #, xticks=([10,50,100],string.([10,50,100])))
yscls = [95e3,140e3]
for (xvals, yvals, yhat) in zip(xvals_4, yvals_4, yscls)
    ysc = yvals[findmin(abs.(ustrip.(xvals) .- yhat))[2]]
    scatter!(p4_varg, ustrip.(xvals) .* 1e-3, yvals ./ ysc , markersize=4, markerstrokewidth=0, xlabel="Distance (km)", ylabel=L"\gamma", label="")
end
# add the fitted model
plot!(p4_varg, [1e-5; ustrip.(xvals_4[1])] .* 1e-3, varg.([1e-5; ustrip.(xvals_4[1])]) ./ varg(150e3), label="Fitted model", lw=3.5, color=:black)
GrIS1980s_DEM.panel_annotate!(p4_varg, "c")
# map
p4_map = plot(aspect_ratio=1, size=(wheight,wwidth), xlims=(-7e5,8e5).*1e-3, ylims=(-3.32e6, -0.78e6).*1e-3, xaxis=false, yaxis=false; palette)
for i_g in is
    scatter!(p4_map, df_not_all.x[i_g].*1e-3, df_not_all.y[i_g].*1e-3, markerstrokewidth=0, markersize=0.5, aspect_ratio=1, label="")
end
GrIS1980s_DEM.panel_annotate!(p4_map, "d")
p_sections = plot(p4_varg, p4_map, layout=grid(1,2,widths=(0.65, 0.35)), size=(wwidth*2, wheight), right_margin=[0mm 10mm], left_margin=[20mm 5mm], bottom_margin=15mm)

plot(p_atm_aero, p_sections, layout=(2,1), size=(wwidth*1.6,wheight*1.8), bottom_margin=[10mm 10mm], top_margin=[20mm 0mm])
savefig(joinpath(fig_dir_main, "fS03.png"))



##############################
# Heatmaps of data, Figure 1 #
##############################

function heatmap_from_df(df_all, sm::Symbol, x, y, dims::Tuple, fname; clims=(-4,4), title)
    # plotting parameters
    Plots.scalefontsizes()
    Plots.scalefontsizes(1.3)
    wsize          = (1000,1200)
    xlabel         = "Easting (km)"
    ylabel         = "Northing (km)"
    colorbar_title = "(-)"
    cmap           = cgrad(:vik, rev=true)
    # plot raster with aerodem data
    id_df_aero = findall(df_all.source .== :aerodem .|| df_all.source .== "aerodem")
    m_plot = zeros(dims)
    m_plot[df_all.idx[id_df_aero]] .= df_all[!,sm][id_df_aero]
    i_nans = m_plot .== 0
    # i_nans[I_no_ocean] .= false
    m_plot[i_nans] .= NaN
    heatmap(x.*1e-3, y.*1e-3, m_plot'; cmap, dpi=300)
    # scatter atm point data
    id_df_atm = findall(df_all.source .== :atm .|| df_all.source .== "atm")
    scatter!(df_all.x[id_df_atm].*1e-3, df_all.y[id_df_atm].*1e-3, marker_z=df_all[!,sm][id_df_atm], label="", markersize=1.2, markerstrokewidth=0, aspect_ratio=1, xlims=(-7e5,8e5).*1e-3, ylims=(-3.32e6, -0.78e6).*1e-3, grid=false, bottommargin=0Plots.cm, rightmargin=(8/5)Plots.cm, leftmargin=(10/5)Plots.cm; wsize, xlabel, ylabel, colorbar_title, title, clims, cmap)
    plot!(outl, fill=nothing, lw=0.2)
    savefig(fname)
    return
end

# plot before standardizing
heatmap_from_df(df_all, :dh, x, y, size(h_ref_m), joinpath(fig_dir_others,"data_non-standardized.png"), clims=(-20,20), title=L"Non-standardized observations $\Delta h$")

# plot after standardizing (FIGURE 2)
heatmap_from_df(df_all, :dh_detrend, x, y, size(h_ref_m), joinpath(fig_dir_main,"f01.png"), clims=(-3.0,3.0), title=L"Standardized observations $\Delta h_\mathrm{std}$")
