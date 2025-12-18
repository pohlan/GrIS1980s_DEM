# define output path for other figures not in paper
fig_dir_others = joinpath("output", "reconstructions", "figures")
mkpath(fig_dir_others)

# download ITS_LIVE velocity
vx_file, vy_file = download_velocity()

# Calculate flowlines from velocity
lstep   = 600.0  # (m), step between two points in flowline
maxdist = 1.2e5  # (m), maximum length of flowline
starting_points = Dict("79N"                     => ( 488378., -1020715.),  # starting coordinates at terminus, determined manually
                       "Daugaard-Jensen"         => ( 559200., -1894248.),
                       "Helheim"                 => ( 313788., -2577999.),
                       "Kangerlussuaq"           => ( 496995., -2295250.),
                       "Kangiata-Nunaata-Sermia" => (-229043., -2818948.),
                       "Petermann"               => (-282347., - 916882.),
                       "Ryder"                   => ( -88571., - 883785.),
                       "Sermeq-Kujalleq"         => (-197941., -2272419.),
                       "Zachariae"               => ( 520829., -1100117.)
)
prof_files = String[]
for (output_name, (px0, py0)) in starting_points
    fname = GrIS1980s_DEM.get_flowline_coords(vx_file, vy_file, lstep, maxdist, px0, py0, output_name)
    push!(prof_files, fname)
end

# helper functions to plot subdomains
function crds_tup(;yb, xl, δx=80000, δy=68000)
    return (;yb,xl,δx,δy)
end
coords = Dict("Helheim"         => crds_tup(yb = -2600000, xl =  260000),
              "Sermeq-Kujalleq" => crds_tup(yb = -2310000, xl = -207000),
              "Humboldt"        => crds_tup(yb = -1080000, xl = -380000),
              "79N"             => crds_tup(yb = -1110715, xl =  405378, δy = 110000, δx = 90000),
              "Zachariae"       => crds_tup(yb = -1180715, xl =  435378, δy = 116000, δx = 95000),
              )
function get_glacier_idx(glacier_name)
    @unpack yb, xl, δx, δy = coords[glacier_name]
    # plot difference map
    ix = findall(xl .< x .< xl + δx)
    iy = findall(yb .< y .< yb + δy)
    return ix, iy
end
rectangle(w, h, x, y) = Shape([x .+ [0,w,w,0], y .+ [0,0,h,h]].*1e-3 ...)

########################################################
# Plot elevations along flowlines, Figures 5a/b and S5 #
########################################################

# files and labels/attributes to zip-loop through
files        = [#get_rec_file_SVD_combined(logλ, r0, grd),
                GrIS1980s_DEM.create_aerodem(grd)[2],
                GrIS1980s_DEM.create_bedmachine_grid(grd)[2],
                # get_rec_file_SVD(logλ, r0, grd),
                # get_rec_file_GP(grd)
                "output/reconstructions/bedmachine1980_SVD_reconstruction_g600.nc",
                "output/reconstructions/bedmachine1980_GP_reconstruction_g600.nc"
                ]
                # "output/reconstructions/rec_kriging_g600_maxn1500.nc" # get_rec_file_kriging(grd, maxn0)
labels       = ["AeroDEM", "GIMP", "SVD method", "GP"]
name_for_col = ["aerodem", "GrIMP", "SVD", "GP"]
bandnm       = ["Band1", "surface", "surface", "surface"]
cols         = Plots.palette(:tol_bright)[2:end]
lstls        = [:solid, :solid, :dot, :dot]
lws          = [6,6,6,6]
z_orders     = [1,1,3,4]
# plotting attributes
Plots.scalefontsizes()
Plots.scalefontsizes(1.9)
ps = Plots.Plot{Plots.GRBackend}[]
xlabel = "Distance along profile (km)"
ylabel = "Surface elevation (m)"
legend = :topleft
# loop through glacier profiles
ixmaxs = Dict()
glacier_names = String[]
for (ip, pf) in enumerate(prof_files)
    # glacier name
    glacier_name = splitext(basename(pf))[1]
    push!(glacier_names, glacier_name)
    # glacier_title = replace(glacier_name, "-" => " ")
    # push!(glacier_titles, glacier_title)
    prof = CSV.read(pf, DataFrame)
    xc = prof.X
    yc = prof.Y
    # initialize figure
    p_i = Plots.plot(size=(GrIS1980s_DEM.wwidth,GrIS1980s_DEM.wheight); xlabel, ylabel, margin = 10Plots.mm,
              fg_legend=:transparent, bg_legend=:transparent, legend)
    # plot the projected elevations of the different DEMs
    for (i,(f, label, col_nm, ls, z_order, band, lw)) in enumerate(zip(files, labels, name_for_col, lstls, z_orders, bandnm, lws))
        dist, vals   = GrIS1980s_DEM.interpolate_raster_to_profile(f, xc, yc; band)
        color = GrIS1980s_DEM.palette_dict[col_nm]
        if label == "AeroDEM"
            i_nonan = findall(.!isnan.(vals))
            i_max   = findmax(vals[i_nonan])[2]
            xmax    = max(min(dist[i_nonan[i_max]]* 1.8, maximum(dist)), 60e3)
            ixmaxs[glacier_name] = findlast(dist .<= xmax)
            xlims=(0, xmax/1e3)
            Plots.plot!(dist./1e3, vals; label, color, ls, lw, z_order, xlims)
        else
            xmax = Plots.xlims(p_i)[2]
            i_nonan = findall(.!isnan.(vals))
            im   = findmin(abs.(xmax .- dist[i_nonan]./1e3))[2]
            ylims=(-1.0, maximum(vals[i_nonan][1:im])+100)
            if label == "Combined SVD/AeroDEM reconstruction"
                dmark = glacier_title == "Helheim" || glacier_title == "Sermeq Kujalleq" ? 3 : 5
                Plots.scatter!(dist[2:dmark:end]./1e3, vals[2:dmark:end]; label, color, ylims, z_order, marker=:diamond, markersize=6.5, markerstrokewidth=0.4, markerstrokecolor="Black", fillalpha=0)
            else
                Plots.plot!(dist./1e3, vals; label, color, ls, lw, ylims, z_order)
            end
            if label == "GP"
                # uncertainty
                dist, vals_std = GrIS1980s_DEM.interpolate_raster_to_profile(f, xc, yc; band = "std_uncertainty")
                vals_std[isnan.(vals_std)] .= 0
                plot!(dist./1e3, vals .- vals_std, fillrange = vals .+ vals_std, label="", fillcolor=color, fillalpha=0.5, lw=0; z_order) #; label, color, ls, lw, ylims, z_order, alpha=0.7)
            end
        end
    end
    push!(ps, p_i)
    Plots.savefig(p_i, joinpath(fig_dir_others, glacier_name*".png"))
end
# # plot all glaciers, Figure S5
# for (i,pp) in enumerate(ps)
#     legend = i == 6 || i == 3
#     ps[i] = plot(pp; legend, markersize=0.3)
# end
# Plots.plot(ps[[2,4,5,6,7,8]]..., size=(2300,2000), left_margin = 12Plots.mm, bottom_margin = 12Plots.mm, topmargin = 12Plots.mm, dpi=300, layout=(3,2))
# Plots.savefig(joinpath(fig_dir_main, "fS05.png"))
# plot two selected glaciers only
name_pos = 0.3
marg_top = 45mm
# Helheim and Sermeq Kujalleq
i_glaciers = findall(glacier_names .== "Helheim" .|| glacier_names .== "Sermeq-Kujalleq")
p1_nolegend = plot(ps[i_glaciers[1]], legend=false, left_margin=30mm, top_margin=marg_top)
GrIS1980s_DEM.panel_annotate!(p1_nolegend, "b") # "  "*glacier_names[i_glaciers[1]]
annotate!(p1_nolegend, mean(xlims(p1_nolegend)), ylims(p1_nolegend)[2]+name_pos*diff([ylims(p1_nolegend)...])[1], text(glacier_names[i_glaciers[1]], 30, :center))
p2          = plot(ps[i_glaciers[2]], legend=:topleft, right_margin=-8mm, left_margin=30mm, top_margin=marg_top)
GrIS1980s_DEM.panel_annotate!(p2, "a")
annotate!(p2, mean(xlims(p2)), ylims(p2)[2]+name_pos*diff([ylims(p2)...])[1], text(glacier_names[i_glaciers[2]], 30, :center))
# outline of Greenland with glacier zoom indicated
p_outl = plot(outl, background_color_inside=nothing, fill=nothing, grid=false, label="", cbar=false, axis=([],false), aspect_ratio=1, lw=0.3)
for glacier in glacier_names[i_glaciers]
    ix, iy = get_glacier_idx(glacier)
    plot!(p_outl, rectangle(x[ix[end]]-x[ix[1]],y[iy[end]]-y[iy[1]],x[ix[1]],y[iy[1]]), fillalpha=0, linewidth=2, linecolor=:red3, label="")
    glacier_title = replace(glacier, "-" => " ")
    xt = (x[ix[end]]-5e4)*1e-3
    yt = (mean(y[iy])+1.5e5)*1e-3
    annotate!([(xt, yt, text(glacier_title, 18))])
end
plot(p_outl, size=(800,900), left_margin=-30mm)
# profiles together with outline
p_profs = plot(p2, p_outl, p1_nolegend, layout=grid(1,3,widths=(0.35,0.2,0.45)), size=(1800,800))

############################################
# Plot difference rec. vs obs, Figure 5c-h #
############################################

# function for plotting zoomed-in map over area
Plots.gr_cbar_width[] = 0.01  # default 0.03
xtick_interval        = 2e4
function plot_map(map_field, glacier_name, panel_letter, flowline_panel, legend_pos::Symbol; cmap=cgrad(:vik, rev=true), clims=(-200,200), id_plot, clabel="   (m)")
    # don't plot anything outside the ice sheet
    map_field[id_plot] .= missing
    # load
    iplot        = 1:ixmaxs[glacier_name]
    prof_name    = prof_files[findfirst(occursin.(glacier_name, prof_files))]
    df           = CSV.read(prof_name, DataFrame)
    @unpack yb, xl, δx, δy = coords[glacier_name]
    # plot map
    ix = findall(xl .< x .< xl + δx)
    iy = findall(yb .< y .< yb + δy)
    xv = Vector(x[ix[1]]:(x[ix[1]]+xtick_interval))
    xtick1 = xv[findfirst(xv .% xtick_interval .== 0)]
    p_dif = heatmap(x[ix].*1e-3, y[iy].*1e-3, map_field[ix,iy]', aspect_ratio=1, xaxis=false, yaxis=false, colorbar_title="", grid=false; cmap, clims)
    annotate!(((xl + 1.32(δx)).*1e-3, (yb + 0.5δy).*1e-3, text(clabel, 18)))
    plot!(p_dif, outl, fill=nothing, xlims=(extrema(x[ix])).*1e-3, ylims=extrema(y[iy]).*1e-3, xticks  = (xtick1:xtick_interval:x[ix[end]]).*1e-3, xtick_direction=:out, lw=0.5)
    # flow line + distance markers + annotation
    plot!(p_dif, df.X[iplot].*1e-3, df.Y[iplot].*1e-3, label="Flowline in ($(flowline_panel))", aspect_ratio=1, lw=3, color=:black, legend=legend_pos, legend_foreground_color=nothing, legend_background_color=:transparent)
    dists, _   = GrIS1980s_DEM.interpolate_raster_to_profile(files[4], df.X, df.Y; band="surface")
    d_marker = [10e3, 30e3, 50e3]
    i_marker = [findmin(abs.(dists .- dm))[2] for dm in d_marker]
    scatter!(p_dif, df.X[i_marker].*1e-3, df.Y[i_marker].*1e-3, color=:slategray, label="", markersize=5, markerstrokewidth=0.5)
    annotate!(p_dif, df.X[i_marker].*1e-3 .+5, df.Y[i_marker].*1e-3 .+4, text.(["10 km", "30 km", "50 km"], 20))
    GrIS1980s_DEM.panel_annotate!(p_dif, panel_letter)
    return p_dif
end
marg_left = 60mm 
marg_bot  = -20mm
title_pos = 26

# difference GP and SVD
h_GP = NCDataset(files[4])["surface"][:,:]
h_SVD = NCDataset(files[3])["surface"][:,:]
mask = NCDataset("data/gris-imbie-1980/outline_mask_g600.nc")["Band1"][:,:]
id_plot = findall(ismissing.(mask))
h_dif = h_GP .- h_SVD
h_aero = NCDataset(files[1])["Band1"][:,:]
# add hashed pattern on top to mark AeroDEM coverage
dd = zeros(size(h_dif)) .+ NaN
dd[.!ismissing.(h_aero)] .= 1
multipolygon = polygonize(==(1.0), range(extrema(x.*1e-3)..., step=1e-3*grd), range(extrema(y.*1e-3)..., step=1e-3*grd), dd)
# Sermeq Kujalleq
p_dif1 = plot_map(h_dif, "Sermeq-Kujalleq", "c", "a", :bottomright; id_plot)
ix, iy = get_glacier_idx("Sermeq-Kujalleq")
plot!(p_dif1, multipolygon, fillstyle=:x, fillcolor=:black, fillalpha=0.3, linewidth=0, label="AeroDEM", right_margin=-35mm, left_margin=marg_left, xlims=xlims(p_dif1))
annotate!(p_dif1, xlims(p_dif1)[1]-title_pos, mean(ylims(p_dif1)), text(L"h_\mathrm{GP}-h_\mathrm{SVD}", "Computer Modern", 30, :center))
# Helheim
p_dif2 = plot_map(h_dif, "Helheim", "d", "b", :bottomleft; id_plot)
ix, iy = get_glacier_idx("Helheim")
plot!(p_dif2, multipolygon, fillstyle=:x, fillcolor=:black, fillalpha=0.3, linewidth=0, label="AeroDEM", margin=10mm, xlims=extrema(x[ix].*1e-3))
# both
p_difs = plot(p_dif1, p_dif2, layout=(1,2), size=(1800,650), bottom_margin=marg_bot)

# GP
p_GP1 = plot_map(h_GP, "Sermeq-Kujalleq", "e", "a", :bottomright, cmap=:terrain, clims=(0,1300); id_plot)
p_GP1 = plot(p_GP1, right_margin=-35mm, left_margin=marg_left)
annotate!(p_GP1, xlims(p_GP1)[1]-title_pos, mean(ylims(p_GP1)), text(L"h_\mathrm{GP}", "Computer Modern", 30, :center))
p_GP2 = plot_map(h_GP, "Helheim", "f", "b", :bottomleft, cmap=:terrain, clims=(0,1900), clabel="         (m)"; id_plot)
p_GP2 = plot(p_GP2, margin=10mm)
p_GPs = plot(p_GP1, p_GP2, layout=(1,2), size=(1800,650), bottom_margin=marg_bot)

# SVD
p_SVD1 = plot_map(h_SVD, "Sermeq-Kujalleq", "g", "a", :bottomright, cmap=:terrain, clims=(0,1300); id_plot)
p_SVD1 = plot(p_SVD1, right_margin=-35mm, left_margin=marg_left)
annotate!(p_SVD1, xlims(p_SVD1)[1]-title_pos, mean(ylims(p_SVD1)), text(L"h_\mathrm{SVD}", "Computer Modern", 30, :center))
p_SVD2 = plot_map(h_SVD, "Helheim", "h", "b", :bottomleft, cmap=:terrain, clims=(0,1900), clabel="         (m)"; id_plot)
p_SVD2 = plot(p_SVD2, margin=10mm)
p_SVDs = plot(p_SVD1, p_SVD2, layout=(1,2), size=(1800,650), bottom_margin=marg_bot)

# all together
plot(p_profs, p_difs, p_GPs, p_SVDs, layout=(4,1), size=(1800,2300))
savefig(joinpath(fig_dir_main, "f05.png"))



############################################
# Plot all-Greenland maps: dif & roughness #
############################################
using Geomorphometry

Plots.scalefontsizes()
Plots.scalefontsizes(1.5)

h_GP = NCDataset(get_rec_file_GP(grd))["surface"][:,:]
h_SVD = NCDataset(get_rec_file_SVD(logλ, r0, grd))["surface"][:,:]
h_dif = h_GP - h_SVD

attr = (;xlims=(-7e5,8e5).*1e-3, ylims=(-3.32e6, -0.78e6).*1e-3, aspect_ratio=1, xaxis=false, yaxis=false, grid=false, left_margin=-5mm, top_margin=0mm, bottom_margin=-10mm)

ps_rough = []
for h in [h_GP, h_SVD]
    rgh = similar(h)
    rgh .= roughness(nomissing(h, 0.0))
    rgh[rgh .==0] .= missing
    p = heatmap(x.*1e-3, y.*1e-3, rgh', clims=(0,25), cmap=cgrad(:lapaz,rev=true), colorbar_title="Roughness"; attr...)
    if isempty(ps_rough)
        panel_annotate_with_title!(p, "a", "   GP")
        # add scalebar and North arrow
        plot!(p, [-7e2,-5e2], [-3e3, -3e3], color=:black, lw=3, label="")
        annotate!(p, [(-6e2, -2.85e3, text("200 km", :center, :Helvetica))])
        plot!(p, [-6e2,-6e2], [-2.6e3, -2.2e3], arrow=true, color=:black, lw=1, label="")
        annotate!(p, [(-6.4e2, -2.45e3, text("N", :right, :Helvetica))])
    else
        GrIS1980s_DEM.panel_annotate!(p, "c")
        panel_annotate_with_title!(p, "c", "   SVD")
    end
    push!(ps_rough, p)
end
p_dif = heatmap(x.*1e-3, y.*1e-3, h_dif', clims=(-50,50), cmap=cgrad(:vik,rev=true), colorbar_title="(m)"; attr...)
panel_annotate_with_title!(p_dif, "b", L"$h_\mathrm{GP}-h_\mathrm{SVD}$")

plot(ps_rough[1], p_dif, ps_rough[2], layout=(1,3), size=(1400,500))
savefig(joinpath(fig_dir_main, "f06.png"))



######################################
# Plot map with flowlines, Figure S6 #
######################################

p = plot(outl, fillalpha=0, aspect_ratio=1, axis=([],false), wsize = (900, 1400), lw=1)
for fname in prof_files
    df            = CSV.read(fname, DataFrame)
    glacier_name  = splitext(basename(fname))[1]
    glacier_title = replace(glacier_name, "-" => "\n")
    iplot         = 1:ixmaxs[glacier_name]
    plot!(df.X[iplot], df.Y[iplot], label="", aspect_ratio=1, lw=3, color="red")
    ann_pos       = df.X[iplot[end]] > 0 ? :right : :left
    dif_x         = df.X[iplot[end]] > 0 ? - 5e3 : 2e4
    dif_y         = -2.5e4
    if glacier_name == "Ryder"
        dif_x = 3e4
        dif_y = 0.0
    end
    Plots.annotate!(df.X[iplot[end]]+dif_x, df.Y[iplot[end]]+dif_y, text(glacier_title, ann_pos, 16, "Computer Modern"))
end
Plots.plot(p, dpi=300)
Plots.savefig(joinpath(fig_dir_main, "fS06.png"))
