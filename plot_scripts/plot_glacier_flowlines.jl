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
files        = [GrIS1980s_DEM.create_aerodem(grd)[2],
                GrIS1980s_DEM.create_bedmachine_grid(grd)[2],
                get_rec_file_SVD(logλ, r0, grd),
                get_rec_file_kriging(grd, maxn0)
                ]
labels       = ["AeroDEM (Korsgaard et al., 2016)", "GrIMP (Howat et al., 2015)", "SVD method", "Kriging"]
name_for_col = ["aerodem", "GrIMP", "SVD", "kriging"]
bandnm       = ["Band1", "surface", "surface", "surface"]
cols         = Plots.palette(:tol_bright)[2:end]
lstls        = [:solid, :solid, :dot, :dot]
lw           = 4
z_orders     = [1,1,1,1]
# plotting attributes
Plots.scalefontsizes()
Plots.scalefontsizes(GrIS1980s_DEM.font_scaling)
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
    p_i = Plots.plot(title="\n"*glacier_title, size=(GrIS1980s_DEM.wwidth,GrIS1980s_DEM.wheight); xlabel, ylabel, margin = 10Plots.mm,
              fg_legend=:transparent, bg_legend=:transparent, legend)
    # plot the projected elevations of the different DEMs
    for (i,(f, label, col_nm, ls, z_order, band)) in enumerate(zip(files, labels, name_for_col, lstls, z_orders, bandnm))
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
            Plots.plot!(dist./1e3, vals; label, color, ls, lw, ylims, z_order)
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
Plots.plot(ps..., size=(2300,2000), left_margin = 12Plots.mm, bottom_margin = 12Plots.mm, topmargin = 12Plots.mm, dpi=300)
Plots.savefig(joinpath(fig_dir_main, "fS05.png"))
# plot two selected glaciers only
i_glaciers = findall(glacier_titles .== "Helheim" .|| glacier_titles .== "Sermeq Kujalleq")
p1_nolegend = plot(ps[i_glaciers[1]], legend=false)
GrIS1980s_DEM.panel_annotate!(p1_nolegend, "b")
p2          = plot(ps[i_glaciers[2]])
GrIS1980s_DEM.panel_annotate!(p2, "a")


############################################
# Plot difference rec. vs obs, Figure 5c/d #
############################################

# load data
dif = nomissing(NCDataset(files[1])[bandnm[1]][:,:] .- NCDataset(files[3])[bandnm[3]][:,:], 0.0)
dif[dif .== 0] .= NaN
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
function plot_dif(glacier_name, panel_letter, flowline_panel, ins_crds::Tuple, arrow_crds::Tuple)
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
    p_dif = heatmap(x[ix], y[iy], dif[ix,iy]', cmap=:RdBu, clims=(-200,200), aspect_ratio=1, size=(1000,900), margin=10Plots.mm, xlabel="Easting (m)", ylabel="Northing (m)", title=" \n \n"*L"$h_\mathrm{obs}- h_\mathrm{SVD}$", colorbar_title="", grid=false)
    annotate!((xl + 1.25δx, yb + 0.5δy, text("(m)", 18)))
    plot!(p_dif, outl, fill=nothing, xlims=extrema(x[ix]), ylims=extrema(y[iy]), xticks  = xtick1:xtick_interval:x[ix[end]], xtick_direction=:out, lw=0.5)
    plot!(p_dif, df.X[iplot], df.Y[iplot], label="Flowline in ($(flowline_panel))", aspect_ratio=1, lw=3, color=:slategray, legend_foreground_color=nothing, legend_background_color=:transparent)
    GrIS1980s_DEM.panel_annotate!(p_dif, panel_letter)
    # insert for Greenland outline
    rectangle(w, h, x, y) = Shape(x .+ [0,w,w,0], y .+ [0,0,h,h])
    ins = bbox(ins_crds[1], ins_crds[2], 0.2, 0.18, :left)
    plot!(p_dif, inset=ins, subplot=2, aspect_ratio=1)
    plot!(p_dif[2], outl, background_color_inside=nothing, fill=nothing, grid=false, label="", cbar=false, axis=([],false), aspect_ratio=1, lw=0.3)
    plot!(p_dif[2], rectangle(x[ix[end]]-x[ix[1]],y[iy[end]]-y[iy[1]],x[ix[1]],y[iy[1]]), fillalpha=0, linewidth=2, linecolor=:red3, label="")
    # draw arrow manually (didn't manage to adjust arrowhead size with arrow=(..))
    x1, xend = x[ix[1]]+arrow_crds[1], x[ix[1]]+arrow_crds[2]
    y1, yend = y[iy[1]]+arrow_crds[3], y[iy[1]]+arrow_crds[4]
    plot!(p_dif[2], [x1,xend], [y1,yend], color=:red3, lw=2, label="")
    plot!(p_dif[2], [xend,xend], [yend,yend+3e5], color=:red3, lw=2, label="")
    plot!(p_dif[2], [xend,xend+sign(arrow_crds[1]-arrow_crds[2])*3e5], [yend,yend+1e5], color=:red3, lw=2, label="")
    return p_dif
end
p_dif1 = plot_dif("Sermeq-Kujalleq", "c", "a", (0.6,0.6), (5e5, 1.5e5, +6e5, +1.5e5))
p_dif2 = plot_dif("Helheim", "d", "b", (0.14,0.6), (-3.5e5, -6e4, +5e5, +1.2e5))
plot(p2, p1_nolegend, p_dif1, p_dif2, wsize=(2500, 1800), left_margin=12Plots.mm, dpi=300)
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
