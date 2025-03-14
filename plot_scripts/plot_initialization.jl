##################################
# Flux divergence maps, Figure 6 #
##################################

# load
sim_folder = joinpath("data", "rec_simulations") # TODO !!
fs = glob(joinpath.(sim_folder, "*.nc"))
grimp_f, svd_f, kriging_f = fs
flux_div_svd     = NCDataset(svd_f)["flux_divergence"][:,:,:]
flux_div_kriging = NCDataset(kriging_f)["flux_divergence"][:,:,:]
flux_0           = findall(flux_div_svd[:,:,1] .== 0)
t                = NCDataset(svd_f)["time"][:]
# coordinates of Sermeq Kujalleq catchment
yb = -2400000
yt = -2140193
xl = -262745
xr =  120000
ix = findall(xl .< x .< xr)
iy = findall(yb .< y .< yt)
# plot
Plots.scalefontsizes()
Plots.scalefontsizes(1.8)
clims = (-30,30)
attr = (; aspect_ratio=1, size=(900,900), colorbar=false, grid=false, clims, cmap=:RdBu, margin=0Plots.mm, link=:all, tick_direction=:out, xticks=[-2e5, -1e5, 0, 1e5])
ps = Plots.Plot{Plots.GRBackend}[]
for (i,flux) in enumerate([flux_div_svd, flux_div_sgs])
    i_t2 = 50
    letter1 = i == 1 ? "a" : "c"
    letter2 = i == 1 ? "b" : "d"
    title1 = i == 1 ? "\n"*L"$t= 1\,\mathrm{a}$" : ""
    year2 = (Year(t[i_t2])-Year(t[1])).value + 1
    title2 = i == 1 ? "\n"*L"$t=%$year2\,\mathrm{a}$" : ""
    xlabel = i == 2 ? "Easting (m)" : ""
    flux[flux_0,:] .= NaN
    # plot at first time step
    p = heatmap(x[ix], y[iy], flux[ix,iy,1]'; title=title1, ylabel="Northing (m)", xlabel="Easting (m)", attr...)
    svd_IceSheetDEM.panel_annotate!(p, letter1)
    plot!(p, outl, xlims=extrema(x[ix]), ylims=extrema(y[iy]), fill=nothing, lw=0.5)
    method = i == 1 ? "SVD method" : "Kriging"
    annotate!(xlims(p)[1] - 3.5e5, mean(ylims(p)), text(method, "Computer Modern", 23, :left))
    if i == 1
        p = plot(p, xtickfontsize=1, xtickfontcolor=:white, xguidefontcolor=:white, bottom_margin=-10Plots.mm, top_margin=20Plots.mm, left_margin=90Plots.mm)
    end
    push!(ps,p)
    # plot later time step
    p = heatmap(x[ix], y[iy], flux[ix,iy,50]'; title=title2, xlabel="Easting (m)", bottom_margin=20Plots.mm, attr...)
    plot!(p, outl, xlims=extrema(x[ix]), ylims=extrema(y[iy]), fill=nothing, lw=0.5)
    if i == 1
        p = plot(p, xtickfontsize=1, xtickfontcolor=:white, xguidefontcolor=:white, bottom_margin=-5Plots.mm)
    end
    p = plot(p, ytickfontsize=1, ytickfontcolor=:white, left_margin=-15Plots.mm)
    svd_IceSheetDEM.panel_annotate!(p, letter2)
    push!(ps,p)
end
# plot insert
rectangle(w, h, x, y) = Shape(x .+ [0,w,w,0], y .+ [0,0,h,h])
ins = bbox(0.65,0.3,0.23,0.42, :left)
plot!(ps[2], inset=ins, subplot=2, aspect_ratio=1)
plot!(ps[2][2], outl, background_color_inside=nothing, fill=nothing, grid=false, label="", cbar=false, axis=([],false))
plot!(ps[2][2], rectangle(x[ix[end]]-x[ix[1]],y[iy[end]]-y[iy[1]],x[ix[1]],y[iy[1]]), fillalpha=0, linewidth=2, linecolor=:darkred, label="")
# plot all panels
p_panels = plot(ps..., layout=(2,2), size=(1800,1000), margin=5Plots.mm)
# colorbar
xx = range(clims...,1000)
zz = zero(xx)' .+ xx
p_c = heatmap(xx, xx, zz, ratio=22, xticks=false, legend=false, fc=cgrad(attr.cmap), lims=clims, framestyle=:box, right_margin=-160Plots.mm, left_margin=-180Plots.mm, top_margin=0Plots.mm, bottom_margin=0Plots.mm, ymirror=true, size=(10,100))
annotate!(300, -15, text("Flux divergence "*L"(\mathrm{m\cdot a}^{-1})", 23, "Computer Modern", :left, :vcenter, rotation=90))
# plot everything together
p_div = plot(p_panels, p_c, layout=(1,2), size=(2000,1000)) #, margin=10Plots.mm)
savefig(joinpath(fig_dir_main, "Figure6.png"))


#########################################
# Mean absolute flux divergence vs time #
#########################################

fig_dir_others = joinpath("output", "reconstructions", "figures")
# mask out non-ice cells
mask = NCDataset(svd_f)["mask"][:,:,1]
id = findall(mask .== 2 .|| mask .== 3)
# mean absolute flux divergence as a function of time
function m_vs_time(A)
    b = zeros(size(A,3))
    for i in axes(A,3)
        b[i] = mean(abs.(A[id,i]))
    end
    return b
end
max_vs_time_svd     = m_vs_time(flux_div_svd)
max_vs_time_kriging = m_vs_time(flux_div_sgs)
# plot
colors = Plots.palette(:batlow10)
p_time = plot(xlabel="Simulation time", ylabel="Mean flux divergence "*L"(\mathrm{a}^{-1})", margin=15Plots.mm)
plot!(t, max_vs_time_svd, label="svd", color=colors[1], lw=3, size=(2000,1000), legend_foreground_color=nothing)
plot!(t, max_vs_time_sgs, label="kriging", lw=3, color=colors[5])
plot!(t, max_vs_time_grimp, label="grimp", lw=3, color=colors[8])
ttck = t[1]:Year(2):t[end]
time_ticks = Dates.format.(ttck,"yyyy")
xticks!(Dates.datetime2epochms.(ttck .- Year(1)), time_ticks)
svd_IceSheetDEM.panel_annotate!(p_time, "a")
savefig(joinpath(fig_dir_others,"mean_flux_div_vs_time.png"))
