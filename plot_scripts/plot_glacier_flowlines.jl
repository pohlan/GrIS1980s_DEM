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
              "Kangerlussuaq"   => crds_tup(yb = -2312000, xl = 401500, δx=110000, δy=93500),
              "Sermeq-Kujalleq" => crds_tup(yb = -2310000, xl = -207000),
              "Ryder"           => crds_tup(yb = -999340, xl = -150571, δx = 140000, δy = 119000),
              "Petermann"        => crds_tup(yb = -1040000, xl = -285000, δx = 160000, δy = 136000),
            #   "Humboldt"        => crds_tup(yb = -1080000, xl = -380000),
              "79N"             => crds_tup(yb = -1110715, xl =  404378, δx = 125000, δy = 106250),
            #   "Zachariae"       => crds_tup(yb = -1180715, xl =  435378, δy = 116000, δx = 95000),
              )
function get_glacier_idx(glacier_name)
    @unpack yb, xl, δx, δy = coords[glacier_name]
    # plot difference map
    ix = findall(xl .< x .< xl + δx)
    iy = findall(yb .< y .< yb + δy)
    return ix, iy
end
rectangle(w, h, x, y) = Shape([x .+ [0,w,w,0], y .+ [0,0,h,h]].*1e-3 ...)

function get_shaded(h, ix, iy; zmin_=NaN, zmax_=NaN)
    zmin = zmax = 0
    h_nomiss = nomissing(h, 0.0)
    hs = hillshade(h_nomiss; cellsize=(grd,grd), azimuth = 315, zenith = 45)
    # normed hillshade
    hs_norm = (hs .- minimum(hs[ix,iy])) ./ (maximum(hs[ix,iy])*0.9 - minimum(hs[ix,iy]))
    hs_norm = hs_norm .^ 0.5
    # normed elevations
    zmin = isnan(zmin_) ? minimum(h_nomiss[ix,iy]) : zmin_
    zmax = isnan(zmax_) ? maximum(h_nomiss[ix,iy]) : zmax_
    znorm = (h_nomiss .- zmin) ./ (zmax*0.9 - zmin)
    znorm = clamp.(znorm, 0, 1)
    # colormap
    cp = cgrad(:terrain)
    colors = get.(Ref(cp), znorm)
    shaded = RGB{Float32}.(
        red.(colors) .* hs_norm,
        green.(colors) .* hs_norm,
        blue.(colors) .* hs_norm
    )
    # mask out values
    shaded[ismissing.(mask)] .= RGB(1,1,1)
    # bla = similar(shaded)
    return shaded, zmin, zmax
end

########################################################
# Plot elevations along flowlines, Figures 6a/b and SXX #
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
z_orders     = [1,1,2,3]
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

name_pos = 0.3
marg_top = 45mm
# function to plot two glaciers
function plot_flowlines_glaciers(glacier1, glacier2)
    i1 = findfirst(glacier_names .== glacier1)
    i2 = findfirst(glacier_names .== glacier2)
    glacier_titles = replace.(glacier_names[[i1,i2]], "-" => " ")
    p1          = plot(ps[i1], legend=:topleft, right_margin=-8mm, left_margin=30mm, top_margin=marg_top)
    GrIS1980s_DEM.panel_annotate!(p1, "a")
    annotate!(p1, mean(xlims(p1)), ylims(p1)[2]+name_pos*diff([ylims(p1)...])[1], text(glacier_titles[1], 30, :center))
    p2_nolegend = plot(ps[i2], legend=false, left_margin=30mm, top_margin=marg_top)
    GrIS1980s_DEM.panel_annotate!(p2_nolegend, "b") # "  "*glacier_names[i_glaciers[1]]
    annotate!(p2_nolegend, mean(xlims(p2_nolegend)), ylims(p2_nolegend)[2]+name_pos*diff([ylims(p2_nolegend)...])[1], text(glacier_titles[2], 30, :center))
    # outline of Greenland with glacier zoom indicated
    p_outl = plot(outl, background_color_inside=nothing, fill=nothing, grid=false, label="", cbar=false, axis=([],false), aspect_ratio=1, lw=0.3)
    for glacier in glacier_names[[i1,i2]]
        ix, iy = get_glacier_idx(glacier)
        plot!(p_outl, rectangle(x[ix[end]]-x[ix[1]],y[iy[end]]-y[iy[1]],x[ix[1]],y[iy[1]]), fillalpha=0, linewidth=2, linecolor=:red3, label="")
        glacier_title = replace(glacier, "-" => " ")
        xt = (x[ix[end]]-5e4)*1e-3
        yt = (mean(y[iy])+1.5e5)*1e-3
        annotate!([(xt, yt, text(glacier_title, 18))])
    end
    plot(p_outl, size=(800,900), left_margin=-30mm)
    # profiles together with outline
    p_profs = plot(p1, p_outl, p2_nolegend, layout=grid(1,3,widths=(0.35,0.2,0.45)), size=(1800,800))
    return p_profs
end

############################################
# Plot difference rec. vs obs, Figure 5c-h #
############################################

# function for plotting zoomed-in map over area
Plots.gr_cbar_width[] = 0.01  # default 0.03
xtick_interval        = 2e4
function plot_map(map_field, glacier_name, panel_letter, flowline_panel, legend_pos::Symbol; cmap=cgrad(:vik, rev=true), clims=nothing, colorbar=true, id_plot, clabel="(m)", no_data=missing, x_shift=5)
    # don't plot anything outside the ice sheet
    map_field[id_plot] .= no_data
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
    yflip = false
    p_dif = heatmap(x[ix].*1e-3, y[iy].*1e-3, map_field[ix,iy]', aspect_ratio=1, xaxis=false, yaxis=false, colorbar_title="", grid=false, yflip=yflip; clims, cmap, colorbar)
    if eltype(map_field) != RGB{Float32}
        annotate!(((xl + 1.32(δx)).*1e-3, (yb + 0.5δy).*1e-3, text(clabel, 18)); yflip)
    end
    plot!(p_dif, outl, fill=nothing, xlims=(extrema(x[ix])).*1e-3, ylims=extrema(y[iy]).*1e-3, xticks  = (xtick1:xtick_interval:x[ix[end]]).*1e-3, xtick_direction=:out, lw=0.5; yflip)
    # flow line + distance markers + annotation
    plot!(p_dif, df.X[iplot].*1e-3, df.Y[iplot].*1e-3, label="Flowline in ($(flowline_panel))", aspect_ratio=1, lw=3, color=:black, legend=legend_pos, legend_foreground_color=nothing, legend_background_color=:transparent; yflip)
    dists, _   = GrIS1980s_DEM.interpolate_raster_to_profile(files[4], df.X, df.Y; band="surface")
    d_marker = 10e3:20e3:(dists[ixmaxs[glacier_name]]-1e3)
    i_marker = [findmin(abs.(dists .- dm))[2] for dm in d_marker]
    scatter!(p_dif, df.X[i_marker].*1e-3, df.Y[i_marker].*1e-3, color=:slategray, label="", markersize=5, markerstrokewidth=0.5; yflip)
    annotate!(p_dif, df.X[i_marker].*1e-3 .+x_shift, df.Y[i_marker].*1e-3 .+4, text.(string.(Int.(d_marker.*1e-3)).* " km", 20); yflip)
    GrIS1980s_DEM.panel_annotate!(p_dif, panel_letter)
    return p_dif
end
marg_left = 60mm
marg_bot  = -5mm

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
function plot_maps_two_glaciers(glacier1, glacier2; title_pos=26, x_shift=5)
    # Sermeq Kujalleq
    p_dif1 = plot_map(h_dif, glacier1, "c", "a", :bottomright, clims=(-200,200); id_plot, x_shift)
    ix, iy = get_glacier_idx(glacier1)
    plot!(p_dif1, multipolygon, fillstyle=:x, fillcolor=:grey, fillalpha=0.4, linewidth=0, label="AeroDEM", right_margin=-25mm, left_margin=marg_left, xlims=xlims(p_dif1))
    annotate!(p_dif1, xlims(p_dif1)[1]-title_pos, mean(ylims(p_dif1)), text(L"h_\mathrm{GP}-h_\mathrm{SVD}", "Computer Modern", 30, :center))
    # add scalebar and North arrow
    plot!(p_dif1, [-145,-135], [-2245, -2245], color=:black, lw=4, label="")
    annotate!(p_dif1, [(-140, -2250, text("10 km", 20, :center))])
    plot!(p_dif1, [-140,-140], [-2265, -2255], color=:black, lw=2, label="")
    plot!(p_dif1, [-138.5,-140], [-2259, -2255], color=:black, lw=2, label="")
    plot!(p_dif1, [-141.5,-140], [-2259, -2255], color=:black, lw=2, label="")
    annotate!(p_dif1, [(-136, -2263, text("N", 22, :bottomright, :Helvetica))])
    # Helheim
    p_dif2 = plot_map(h_dif, glacier2, "d", "b", :bottomleft, clims=(-200,200); id_plot, x_shift)
    ix, iy = get_glacier_idx(glacier2)
    plot!(p_dif2, multipolygon, fillstyle=:x, fillcolor=:grey, fillalpha=0.4, linewidth=0, label="AeroDEM", margin=10mm, xlims=extrema(x[ix].*1e-3))
    # add scalebar and North arrow
    plot!(p_dif2, [327,337], [-2598, -2598], color=:black, lw=4, label="")
    annotate!(p_dif2, [(332, -2593, text("10 km", 20, :center))])
    plot!(p_dif2, [332,332], [-2589, -2579], color=:black, lw=2, label="")
    plot!(p_dif2, [330.5,332], [-2583, -2579], color=:black, lw=2, label="")
    plot!(p_dif2, [333.5,332], [-2583, -2579], color=:black, lw=2, label="")
    annotate!(p_dif2, [(336, -2587, text("N", 22, :bottomleft, :Helvetica))])
    # both
    p_difs = plot(p_dif1, p_dif2, layout=(1,2), size=(1800,650), bottom_margin=marg_bot)

    # GP
    shaded, zmin_, zmax_ = get_shaded(h_SVD, get_glacier_idx(glacier1)...)
    shaded, _ = get_shaded(h_GP, get_glacier_idx(glacier1)...; zmin_, zmax_)
    p_GP1 = plot_map(shaded, glacier1, "e", "a", :bottomright, colorbar=false, cmap=nothing, no_data=RGB(1,1,1); id_plot, x_shift)
    p_GP1 = plot(p_GP1, right_margin=-35mm, left_margin=marg_left-110mm)
    annotate!(p_GP1, xlims(p_GP1)[1]-title_pos, mean(ylims(p_GP1)), text(L"h_\mathrm{GP}", "Computer Modern", 30, :center))
    shaded, zmin_, zmax_ = get_shaded(h_SVD, get_glacier_idx(glacier2)...)
    shaded, _ = get_shaded(h_GP, get_glacier_idx(glacier2)...; zmin_, zmax_)
    p_GP2 = plot_map(shaded, glacier2, "f", "b", :bottomleft, colorbar=false, cmap=nothing, no_data=RGB(1,1,1); id_plot, x_shift)
    p_GP2 = plot(p_GP2, margin=10mm, left_margin=-95mm)
    p_GPs = plot(p_GP1, p_GP2, layout=(1,2), size=(1800,650), bottom_margin=marg_bot)

    # SVD
    shaded, _ = get_shaded(h_SVD, get_glacier_idx(glacier1)...)
    p_SVD1 = plot_map(shaded, glacier1, "g", "a", :bottomright, colorbar=false, cmap=nothing, no_data=RGB(1,1,1); id_plot, x_shift)
    p_SVD1 = plot(p_SVD1, right_margin=-35mm, left_margin=marg_left-110mm)
    annotate!(p_SVD1, xlims(p_SVD1)[1]-title_pos, mean(ylims(p_SVD1)), text(L"h_\mathrm{SVD}", "Computer Modern", 30, :center))
    shaded, _ = get_shaded(h_SVD, get_glacier_idx(glacier2)...)
    p_SVD2 = plot_map(shaded, glacier2, "h", "b", :bottomleft, colorbar=false, cmap=nothing, no_data=RGB(1,1,1); id_plot, x_shift)
    p_SVD2 = plot(p_SVD2, margin=10mm, left_margin=-95mm)
    p_SVDs = plot(p_SVD1, p_SVD2, layout=(1,2), size=(1800,650), bottom_margin=marg_bot)
    return p_difs, p_GPs, p_SVDs
end


# all together for Helheim and Sermeq Kujalleq
p_profs = plot_flowlines_glaciers("Sermeq-Kujalleq", "Helheim")
p_difs, p_GPs, p_SVDs = plot_maps_two_glaciers("Sermeq-Kujalleq", "Helheim")
plot(p_profs, p_difs, p_GPs, p_SVDs, layout=(4,1), size=(1800,2300))
savefig(joinpath(fig_dir_main, "f05.png"))

# 79N and Kangerlussuaq
p_profs = plot_flowlines_glaciers("79N", "Ryder")
p_difs, p_GPs, p_SVDs = plot_maps_two_glaciers("79N", "Ryder", title_pos=40, x_shift=17)
plot(p_profs, p_difs, p_GPs, p_SVDs, layout=(4,1), size=(1800,2300))
savefig(joinpath(fig_dir_main, "fSXX_79N_Ryder.png"))

# Ryder and Kangerlussuaq
p_profs = plot_flowlines_glaciers("Petermann", "Kangerlussuaq")
p_difs, p_GPs, p_SVDs = plot_maps_two_glaciers("Petermann", "Kangerlussuaq", title_pos=40, x_shift=17)
plot(p_profs, p_difs, p_GPs, p_SVDs, layout=(4,1), size=(1800,2300))
savefig(joinpath(fig_dir_main, "fSXX_Petermann_Kangerlussuaq.png"))


######################################################
# Plot all-Greenland maps: dif & roughness, Figure 7 #
######################################################

Plots.scalefontsizes()
Plots.scalefontsizes(1.5)

h_GP = NCDataset(get_rec_file_GP(grd))["surface"][:,:]
h_SVD = NCDataset(get_rec_file_SVD(logλ, r0, grd))["surface"][:,:]
h_dif = h_GP - h_SVD

attr = (;xlims=(-7e5,8e5).*1e-3, ylims=(-3.32e6, -0.78e6).*1e-3, aspect_ratio=1, xaxis=false, yaxis=false, grid=false, left_margin=-5mm, top_margin=0mm, bottom_margin=-10mm)

ps_rough = []
for h in [h_GP, h_SVD]
    rgh = similar(h)
    rgh .= TRI(nomissing(h, 0.0))
    rgh[rgh .==0] .= missing
    p = heatmap(x.*1e-3, y.*1e-3, rgh', clims=(0,45), cmap=cgrad(:lapaz,rev=true), colorbar_title="TRI (m)"; attr...)
    if isempty(ps_rough)
        GrIS1980s_DEM.panel_annotate_with_title!(p, "a", "   GP")
        # add scalebar and North arrow
        plot!(p, [-7e2,-5e2], [-3e3, -3e3], color=:black, lw=3, label="")
        annotate!(p, [(-6e2, -2.85e3, text("200 km", :center, :Helvetica))])
        plot!(p, [-6e2,-6e2], [-2.6e3, -2.2e3], arrow=true, color=:black, lw=1, label="")
        annotate!(p, [(-6.4e2, -2.45e3, text("N", :right, :Helvetica))])
    else
        GrIS1980s_DEM.panel_annotate!(p, "c")
        GrIS1980s_DEM.panel_annotate_with_title!(p, "c", "   SVD")
    end
    push!(ps_rough, p)
end
p_dif = heatmap(x.*1e-3, y.*1e-3, h_dif', clims=(-50,50), cmap=cgrad(:vik,rev=true), colorbar_title="(m)"; attr...)
GrIS1980s_DEM.panel_annotate_with_title!(p_dif, "b", L"$h_\mathrm{GP}-h_\mathrm{SVD}$")

plot(ps_rough[1], p_dif, ps_rough[2], layout=(1,3), size=(1400,500))
savefig(joinpath(fig_dir_main, "f07.png"))
