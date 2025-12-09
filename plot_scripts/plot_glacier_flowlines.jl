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

########################################################
# Plot elevations along flowlines, Figures 5a/b and S5 #
########################################################

# files and labels/attributes to zip-loop through
files        = [get_rec_file_SVD_combined(logλ, r0, grd),
                GrIS1980s_DEM.create_aerodem(grd)[2],
                GrIS1980s_DEM.create_bedmachine_grid(grd)[2],
                get_rec_file_SVD(logλ, r0, grd),
                get_rec_file_GP(grd)]
                # "output/reconstructions/rec_kriging_g600_maxn1500.nc" # get_rec_file_kriging(grd, maxn0)

labels       = ["Combined SVD/AeroDEM reconstruction", "AeroDEM (Korsgaard et al., 2016)", "GrIMP (Howat et al., 2015)", "SVD method", "GP"]
name_for_col = ["combined", "aerodem", "GrIMP", "SVD", "GP"]
bandnm       = ["surface", "Band1", "surface", "surface", "surface"]
cols         = Plots.palette(:tol_bright)[2:end]
lstls        = [:dashdotdot, :solid, :solid, :dot, :dot]
lws          = [3,4,4,4,4]
z_orders     = [1,1,1,1,1]
# plotting attributes
Plots.scalefontsizes()
Plots.scalefontsizes(1.9)
ps = Plots.Plot{Plots.GRBackend}[]
xlabel = "Distance along profile (km)"
ylabel = "Surface elevation (m)"
legend = :topleft
# loop through glacier profiles
ixmaxs = Dict()
glacier_titles = String[]
for (ip, pf) in enumerate(prof_files)
    # glacier name
    glacier_name = splitext(basename(pf))[1]
    glacier_title = replace(glacier_name, "-" => " ")
    push!(glacier_titles, glacier_title)
    prof = CSV.read(pf, DataFrame)
    xc = prof.X
    yc = prof.Y
    # initialize figure
    p_i = Plots.plot(title=glacier_title, size=(GrIS1980s_DEM.wwidth,GrIS1980s_DEM.wheight); xlabel, ylabel, margin = 10Plots.mm,
              fg_legend=:transparent, bg_legend=:transparent, legend)
    # plot the projected elevations of the different DEMs
    for (i,(f, label, col_nm, ls, z_order, band, lw)) in enumerate(zip(files, labels, name_for_col, lstls, z_orders, bandnm, lws))
        dist, vals   = GrIS1980s_DEM.interpolate_raster_to_profile(f, xc, yc; band)
        color = GrIS1980s_DEM.palette_dict[col_nm]
        if label == "AeroDEM (Korsgaard et al., 2016)"
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
            if label == "Kriging" || label == "SVD method"
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
# plot all glaciers, Figure S5
for (i,pp) in enumerate(ps)
    legend = i == 6 || i == 3
    ps[i] = plot(pp; legend, markersize=0.3)
end
Plots.plot(ps[[2,4,5,6,7,8]]..., size=(2300,2000), left_margin = 12Plots.mm, bottom_margin = 12Plots.mm, topmargin = 12Plots.mm, dpi=300, layout=(3,2))
Plots.savefig(joinpath(fig_dir_main, "fS05.png"))
# plot two selected glaciers only
i_glaciers = findall(glacier_titles .== "Helheim" .|| glacier_titles .== "Sermeq Kujalleq")
p1_nolegend = plot(ps[i_glaciers[1]], legend=false, aspect_ratio=0.02)
GrIS1980s_DEM.panel_annotate!(p1_nolegend, "b")
p2          = plot(ps[i_glaciers[2]], aspect_ratio=0.02, legend=:topleft)
GrIS1980s_DEM.panel_annotate!(p2, "a")


############################################
# Plot difference rec. vs obs, Figure 5c/d #
############################################

# load data
# dif = nomissing(NCDataset(files[2])[bandnm[2]][:,:] .- NCDataset(files[4])[bandnm[4]][:,:], 0.0)
# dif[dif .== 0] .= NaN
# coordinates of glacier catchments
function crds_tup(;yb, xl, δx=80000, δy=68000)
    return (;yb,xl,δx,δy)
end
coords = Dict("Helheim"         => crds_tup(yb = -2600000, xl =  260000),
              "Sermeq-Kujalleq" => crds_tup(yb = -2310000, xl = -227000),
              "Humboldt"        => crds_tup(yb = -1080000, xl = -380000),
              "79N"             => crds_tup(yb = -1110715, xl =  405378, δy = 110000, δx = 90000),
              "Zachariae"       => crds_tup(yb = -1180715, xl =  435378, δy = 116000, δx = 95000),
              )
# plot
Plots.gr_cbar_width[] = 0.01  # default 0.03
xtick_interval        = 2e4
function plot_dif(map_field, glacier_name, panel_letter, flowline_panel, legend_pos::Symbol, shift_x=0; title, cmap=cgrad(:vik, rev=true), clims=(-200,200), inset=false, ins_crds=(0,0), arrow_crds=(0,0,0,0))
    # load
    iplot        = 1:ixmaxs[glacier_name]
    prof_name    = prof_files[findfirst(occursin.(glacier_name, prof_files))]
    df           = CSV.read(prof_name, DataFrame)
    @unpack yb, xl, δx, δy = coords[glacier_name]
    # plot difference map
    ix = findall(xl .< x .< xl + δx)
    iy = findall(yb .< y .< yb + δy)
    xv = Vector(x[ix[1]]:(x[ix[1]]+xtick_interval))
    xtick1 = xv[findfirst(xv .% xtick_interval .== 0)]
    p_dif = heatmap(x[ix].*1e-3, y[iy].*1e-3, map_field[ix,iy]', aspect_ratio=1, size=(1000,900), xlabel="Easting (km)", ylabel="Northing (km)", colorbar_title="", grid=false; cmap, clims, title)
    annotate!(((xl + 1.32(δx+shift_x)).*1e-3, (yb + 0.5δy).*1e-3, text("(m)", 18)))
    plot!(p_dif, outl, fill=nothing, xlims=(extrema(x[ix]).+shift_x).*1e-3, ylims=extrema(y[iy]).*1e-3, xticks  = (xtick1:xtick_interval:x[ix[end]]).*1e-3, xtick_direction=:out, lw=0.5)

    # flow line + distance markers + annotation
    plot!(p_dif, df.X[iplot].*1e-3, df.Y[iplot].*1e-3, label="Flowline in ($(flowline_panel))", aspect_ratio=1, lw=3, color=:black, legend=legend_pos, legend_foreground_color=nothing, legend_background_color=:transparent)
    dists, _   = GrIS1980s_DEM.interpolate_raster_to_profile(files[4], df.X, df.Y; band="surface")
    d_marker = [10e3, 30e3, 50e3]
    i_marker = [findmin(abs.(dists .- dm))[2] for dm in d_marker]
    scatter!(p_dif, df.X[i_marker].*1e-3, df.Y[i_marker].*1e-3, color=:slategray, label="", markersize=5, markerstrokewidth=0.5)
    annotate!(p_dif, df.X[i_marker].*1e-3 .+2, df.Y[i_marker].*1e-3 .+3, text.(["10 km", "30 km", "50 km"], 20))
    GrIS1980s_DEM.panel_annotate!(p_dif, panel_letter)
    if inset
        # insert for Greenland outline
        rectangle(w, h, x, y) = Shape([x .+ [0,w,w,0], y .+ [0,0,h,h]].*1e-3 ...)
        ins = bbox(ins_crds[1], ins_crds[2], 0.2, 0.4, :left)
        plot!(p_dif, inset=ins, subplot=2, aspect_ratio=1)
        plot!(p_dif[2], outl, background_color_inside=nothing, fill=nothing, grid=false, label="", cbar=false, axis=([],false), aspect_ratio=1, lw=0.3)
        plot!(p_dif[2], rectangle(x[ix[end]]-x[ix[1]],y[iy[end]]-y[iy[1]],x[ix[1]],y[iy[1]]), fillalpha=0, linewidth=2, linecolor=:red3, label="")
        # draw arrow manually (didn't manage to adjust arrowhead size with arrow=(..))
        x1, xend = (x[ix[1]]+arrow_crds[1], x[ix[1]]+arrow_crds[2]).*1e-3
        y1, yend = (y[iy[1]]+arrow_crds[3], y[iy[1]]+arrow_crds[4]).*1e-3
        plot!(p_dif[2], [x1,xend], [y1,yend], color=:red3, lw=2, label="")
        plot!(p_dif[2], [xend,xend], [yend,yend+3e2], color=:red3, lw=2, label="")
        plot!(p_dif[2], [xend,xend+sign(arrow_crds[1]-arrow_crds[2])*3e2], [yend,yend+1e2], color=:red3, lw=2, label="")
    end
    return p_dif
end

# GP
h_GP = NCDataset(files[5])["surface"][:,:]
p_GP1 = plot_dif(h_GP, "Sermeq-Kujalleq", "e", "a", :bottomright, 2e3, title=" \n"*L"$h_\mathrm{GP}$", cmap=:terrain, clims=(0,1300))
p_GP2 = plot_dif(h_GP, "Helheim", "f", "b", :bottomright, -1.5e3, title=" \n"*L"$h_\mathrm{GP}$", cmap=:terrain, clims=(0,1900))

# SVD
h_SVD = NCDataset(files[4])["surface"][:,:]
p_SVD1 = plot_dif(h_SVD, "Sermeq-Kujalleq", "g", "a", :bottomright, 2e3, title=" \n"*L"$h_\mathrm{SVD}$", cmap=:terrain, clims=(0,1300))
p_SVD2 = plot_dif(h_SVD, "Helheim", "h", "b", :bottomright, -1.5e3, title=" \n"*L"$h_\mathrm{SVD}$", cmap=:terrain, clims=(0,1900))

# difference GP and SVD
h_dif = h_GP .- h_SVD
h_aero = NCDataset(files[2])["Band1"][:,:]
p_dif1 = plot_dif(h_dif, "Sermeq-Kujalleq", "g", "a", :bottomright, 2e3, title=" \n"*L"$h_\mathrm{SVD}$")
dd = zeros(size(h_dif)) .+ NaN
dd[.!ismissing.(h_aero)] .= 1
xl = xlims(p_dif1); ix = findall(xl[1] .<= x.*1e-3 .<= xl[2])
yl = ylims(p_dif1); iy = findall(yl[1] .<= y.*1e-3 .<= yl[2])
heatmap!(p_dif1, x[ix].*1e-3, y[iy].*1e-3, dd[ix,iy]', alpha=0.2, color=:gray, colorbar=false, xlims=xl, ylims=yl)
p_dif2 = plot_dif(h_dif, "Helheim", "h", "b", :bottomright, -1.5e3, title=" \n"*L"$h_\mathrm{SVD}$")
xl = xlims(p_dif2); ix = findall(xl[1] .<= x.*1e-3 .<= xl[2])
yl = ylims(p_dif2); iy = findall(yl[1] .<= y.*1e-3 .<= yl[2])
heatmap!(p_dif2, x[ix].*1e-3, y[iy].*1e-3, dd[ix,iy]', alpha=0.2, color=:gray, colorbar=false, xlims=xl, ylims=yl)


plot(p2, p1_nolegend, p_dif1, p_dif2, p_GP1, p_GP2, p_SVD1, p_SVD2, wsize=(2400, 2900), dpi=300, layout=grid(4, 2, widths=(0.5, 0.5), heights=(0.16, 0.28,0.28,0.28)), left_margin=15mm)
savefig(joinpath(fig_dir_main, "f05.png"))

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
