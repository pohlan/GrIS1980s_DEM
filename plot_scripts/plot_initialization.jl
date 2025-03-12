using svd_IceSheetDEM, NCDatasets, Plots, StatsBase, Glob, Shapefile, Dates, LaTeXStrings

outline_shp_file = joinpath("data", "gris-imbie-1980", "gris-outline-imbie-1980_updated_crs.shp")

sim_folder = joinpath("data", "rec_simulations")
fs = glob(joinpath.(sim_folder, "*.nc"))
grimp_f, svd_f, sgs_f = fs

ds_sgs = NCDataset(sgs_f)
ds_grimp = NCDataset(grimp_f)
ds_svd = NCDataset(svd_f)
x      = ds_svd["x"][:]
y      = ds_svd["y"][:]
t      = ds_svd["time"][:]

shp              = Shapefile.shapes(Shapefile.Table(outline_shp_file))

Plots.scalefontsizes()
Plots.scalefontsizes(1.8)

# coordinates of Sermeq Kujalleq catchment
yb = -2400000
yt = -2140193
xl = -262745
xr =  120000

ix = findall(xl .< x .< xr)
iy = findall(yb .< y .< yt)

nm = ["flux_divergence", "dHdt", "usurf"][1]
flux_div_svd = ds_svd[nm][:,:,:]
flux_div_grimp = ds_grimp[nm][:,:,:]
flux_div_sgs = ds_sgs[nm][:,:,:]

clims = (-30,30)
attr = (; aspect_ratio=1, size=(900,900), colorbar=false, grid=false, clims, cmap=:RdBu, margin=0Plots.mm, link=:all, tick_direction=:out, xticks=[-2e5, -1e5, 0, 1e5])

flux_0 = findall(flux_div_svd[:,:,1] .== 0)

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
    p = heatmap(x[ix], y[iy], flux[ix,iy,1]'; title=title1, ylabel="Northing (m)", xlabel="Easting (m)", attr...)
    svd_IceSheetDEM.panel_annotate!(p, letter1)
    plot!(p, shp, xlims=extrema(x[ix]), ylims=extrema(y[iy]), fill=nothing, lw=0.5)
    method = i == 1 ? "SVD method" : "Kriging"
    annotate!(xlims(p)[1] - 3.5e5, mean(ylims(p)), text(method, "Computer Modern", 23, :left))
    if i == 1
        p = plot(p, xtickfontsize=1, xtickfontcolor=:white, xguidefontcolor=:white, bottom_margin=-10Plots.mm, top_margin=20Plots.mm, left_margin=90Plots.mm) #
    end
    push!(ps,p)
    p = heatmap(x[ix], y[iy], flux[ix,iy,50]'; title=title2, xlabel="Easting (m)", bottom_margin=20Plots.mm, attr...)
    plot!(p, shp, xlims=extrema(x[ix]), ylims=extrema(y[iy]), fill=nothing, lw=0.5)
    if i == 1
        p = plot(p, xtickfontsize=1, xtickfontcolor=:white, xguidefontcolor=:white, bottom_margin=-5Plots.mm)
    end
    p = plot(p, ytickfontsize=1, ytickfontcolor=:white, left_margin=-15Plots.mm)
    svd_IceSheetDEM.panel_annotate!(p, letter2)
    push!(ps,p)
end

# plot as insert to previous plot p_map
rectangle(w, h, x, y) = Shape(x .+ [0,w,w,0], y .+ [0,0,h,h])
ins = bbox(0.65,0.3,0.23,0.42, :left)
plot!(ps[2], inset=ins, subplot=2, aspect_ratio=1)
plot!(ps[2][2], shp, background_color_inside=nothing, fill=nothing, grid=false, label="", cbar=false, axis=([],false))
plot!(ps[2][2], rectangle(x[ix[end]]-x[ix[1]],y[iy[end]]-y[iy[1]],x[ix[1]],y[iy[1]]), fillalpha=0, linewidth=2, linecolor=:darkred, label="")

# plot all panels
p_panels = plot(ps..., layout=(2,2), size=(1800,1000), margin=5Plots.mm)
savefig("bla2.png")

xx = range(clims...,1000)
zz = zero(xx)' .+ xx
p_c = heatmap(xx, xx, zz, ratio=22, xticks=false, legend=false, fc=cgrad(attr.cmap), lims=clims, framestyle=:box, right_margin=-160Plots.mm, left_margin=-180Plots.mm, top_margin=0Plots.mm, bottom_margin=0Plots.mm, ymirror=true, size=(10,100))
annotate!(300, -15, text("Flux divergence "*L"(\mathrm{m\cdot a}^{-1})", 23, "Computer Modern", :left, :vcenter, rotation=90))
plot(p_c)

# plot again everything together
# layout = @layout [a{0.9w,1.0h} b{0.7h}]

p_div = plot(p_panels, p_c, layout=(1,2), size=(2000,1000)) #, margin=10Plots.mm)
savefig("Flux_divergence_comparison.png")

# p1 = heatmap(x[ix], y[iy], flux_div_svd[ix,iy,tsp]',   title="SVD";     attr...)
# p2 = heatmap(x[ix], y[iy], flux_div_sgs[ix,iy,tsp]',   title="kriging"; attr...)
# p3 = heatmap(x[ix], y[iy], flux_div_grimp[ix,iy,tsp]', title="grimp";   attr...)
# p4 = plot(1:4)
# plot(p1,p2,p3,p4, layout=(3,2), size=(900, 900))

# tsp = 131
# clm = 20
# p1 = heatmap(x[ix], y[iy], flux_div_svd[ix,iy,tsp]', cmap=:RdBu, clims=(-clm,clm),   title="SVD", aspect_ratio=1, size=(800,900))
# p2 = heatmap(x[ix], y[iy], flux_div_sgs[ix,iy,tsp]', cmap=:RdBu, clims=(-clm,clm),   title="kriging", aspect_ratio=1, size=(800,900))
# p3 = heatmap(x[ix], y[iy], flux_div_grimp[ix,iy,tsp]', cmap=:RdBu, clims=(-clm,clm), title="grimp", margin=5Plots.mm, aspect_ratio=1, size=(900,900))
# plot(p1,p2,p3, layout=(1,3), size=(1800, 800))
# savefig("bla2.png")



################
# Mean vs time #
################

mask = ds_sgs["mask"][:,:,1]
id = findall(mask .== 2 .|| mask .== 3)

# plot as a function of time
function m_vs_time(A, f, ds)
    b = zeros(size(A,3))
    for i in axes(A,3)
        # mask = ds["mask"][:,:,i]
        b[i] = f(abs.(A[id,i]))
    end
    return b
end
max_vs_time_svd = m_vs_time(flux_div_svd,     mean, ds_svd)
max_vs_time_sgs = m_vs_time(flux_div_sgs,     mean, ds_sgs)
max_vs_time_grimp = m_vs_time(flux_div_grimp, mean, ds_grimp)

colors = Plots.palette(:batlow10)

p_time = plot(xlabel="Simulation time", ylabel="Mean flux divergence "*L"(\mathrm{a}^{-1})", margin=15Plots.mm)
plot!(t, max_vs_time_svd, label="svd", color=colors[1], lw=3, size=(2000,1000), legend_foreground_color=nothing)
plot!(t, max_vs_time_sgs, label="kriging", lw=3, color=colors[5])
plot!(t, max_vs_time_grimp, label="grimp", lw=3, color=colors[8])
ttck = t[1]:Year(2):t[end]
time_ticks = Dates.format.(ttck,"yyyy")
xticks!(Dates.datetime2epochms.(ttck .- Year(1)), time_ticks)
svd_IceSheetDEM.panel_annotate!(p_time, "a")
savefig("bla.png")

plot(p_time, p_div, layout=(1,2), size=(2000, 800))
savefig("bla2.png")



# Plots.scalefontsizes()
# Plots.scalefontsizes(1.6)
# @unpack maxns, m_interps = load("output/validation/kriging_findmaxn.jld2")
# ps = Plots.Plot{Plots.GRBackend}[]
# clims=(-4,4)
# cmap = :RdBu
# for m in axes(m_interps, 3)
#     # cbar = m == 3 ? true : false
#     yticks = m == 1 ? true : false
#     pi = heatmap(x[xsp], y[ysp], m_interps[:,:,m]', title="\n"*L"\mathrm{n_{obs}} = "*"$(maxns[m])", aspect_ratio=1, wsize=(1500,400), grid=false, cbar=false, tick_direction=:out, titlefontsize=18; clims, cmap)
#     if m !== 1
#         pi = plot(pi, ytickfontsize=1, ytickfontcolor=:white)
#     end
#     plot!(shp, xlims=extrema(x[xsp]), ylims=extrema(y[ysp]), fill=nothing, lw=0.5)
#     push!(ps, pi)
# end
# p_panels = plot(ps..., size=(3000, 500), layout=(1,4), leftmargin=10Plots.mm, rightmargin=10Plots.mm, topmargin=-10Plots.mm, bottommargin=-10Plots.mm)
