using svd_IceSheetDEM, NCDatasets, Glob, CSV, DataFrames, Shapefile, Plots, UnPack, LaTeXStrings

#####################################
# calculate flowlines from velocity #
#####################################

vx_file = "data/ITS_LIVE_velocity_120m_RGI05A_0000_v02_vx.nc"
vy_file = "data/ITS_LIVE_velocity_120m_RGI05A_0000_v02_vy.nc"
lstep   = 600.0  # (m), step between two points in flowline
maxdist = 1.2e5  # (m), maximum length of flowline

starting_points = Dict("79N"                     => ( 488378., -1020715.),
                       "Daugaard-Jensen"         => ( 559200., -1894248.),
                       "Helheim"                 => ( 313788., -2577999.),
                 #    "Humboldt"                => (),
                       "Kangerlussuaq"           => ( 496995., -2295250.),
                       "Kangiata-Nunaata-Sermia" => (-229043., -2818948.),
                       "Petermann"               => (-282347., - 916882.),
                       "Ryder"                   => ( -88571., - 883785.),
                       "Sermeq-Kujalleq"         => (-197941., -2272419.),
                    #    "Storstrommen"            => (),
                       "Zachariae"               => ( 520829., -1100117.)
)

prof_files = String[]
for (output_name, (px0, py0)) in starting_points
    fname = svd_IceSheetDEM.get_flowline_coords(vx_file, vy_file, lstep, maxdist, px0, py0, output_name)
    push!(prof_files, fname)
end

#####################################
# plot elevations along flowlines   #
#####################################

fig_dir = joinpath("output", "profiles", "figures")
mkpath(fig_dir)

# randn_files = glob("output/rec_rand_id_*.nc")
# SEQ_files    = glob("output/SEQ/simulations/SEQ_maxngh_100_g600_id_*_upsampled.nc")
# SVD_RF_files = glob("output/SVD_RF/simulations/rec_rand_id_*.nc")
files        = ["data/aerodem/aerodem_rm-filtered_geoid-corr_g600.nc",
                "data/bedmachine/bedmachine_g600.nc",
                "output/SVD_reconstruction/rec_lambda_1e7_g600_r500_h.nc",
                "output/geostats_interpolation/kriging/rec_kriging_g600_maxn1600.nc"
               ]
labels = ["Korsgaard et al., 2016", "GrIMP (Howat et al., 2015)", "SVD method", "Kriging"]

name_for_col =  ["aerodem", "GrIMP", "SVD_h", "kriging"]
bandnm = ["Band1", "surface", "surface", "surface"]
cols   = Plots.palette(:tol_bright)[2:end]
lstls  = [:solid, :solid, :dot, :dot]
lw     = 4
z_orders = [1,1,1,1]

Plots.scalefontsizes()
Plots.scalefontsizes(svd_IceSheetDEM.font_scaling)
ps = Plots.Plot{Plots.GRBackend}[]

ixmaxs = Dict()
i20    = Dict()
glacier_titles = String[]
for (ip, pf) in enumerate(prof_files) #[[3,4,5,6]])
    glacier_name = splitext(basename(pf))[1]
    glacier_title = replace(glacier_name, "-" => " ")
    push!(glacier_titles, glacier_title)
    prof = CSV.read(pf, DataFrame)
    xc = prof.X
    yc = prof.Y

    # set up figure
    xlabel = "Distance along profile (km)"
    ylabel = "Surface elevation (m)"
    legend = :topleft
    p_i = Plots.plot(title="\n"*glacier_title, size=(svd_IceSheetDEM.wwidth,svd_IceSheetDEM.wheight); xlabel, ylabel, margin = 10Plots.mm,
              fg_legend=:transparent, bg_legend=:transparent, legend)

    # plot GrIMP, reconstruction and aerodem
    for (i,(f, label, col_nm, ls, z_order, band)) in enumerate(zip(files, labels, name_for_col, lstls, z_orders, bandnm))
        dist, vals   = svd_IceSheetDEM.interpolate_raster_to_profile(f, xc, yc; band)
        color = svd_IceSheetDEM.palette_dict[col_nm]
        if label == "Korsgaard et al., 2016"
            i_nonan = findall(.!isnan.(vals))
            i_max   = findmax(vals[i_nonan])[2]
            xmax    = max(min(dist[i_nonan[i_max]]* 1.8, maximum(dist)), 60e3)
            ixmaxs[glacier_name] = findlast(dist .<= xmax)
            i20[glacier_name] = findlast(dist .<= 15e3)
            xlims=(0, xmax/1e3)
            Plots.plot!(dist./1e3, vals; label, color, ls, lw, z_order, xlims)
            # Plots.scatter!([dist[i20[glacier_name]]./1e3], [vals[i20[glacier_name]]], color="black", markersize=6, label="")
        else
            xmax = Plots.xlims(p_i)[2]
            i_nonan = findall(.!isnan.(vals))
            im   = findmin(abs.(xmax .- dist[i_nonan]./1e3))[2]
            ylims=(-1.0, maximum(vals[i_nonan][1:im])+100)
            Plots.plot!(dist./1e3, vals; label, color, ls, lw, ylims, z_order)
            if label == "Kriging" || label == "SVD method"
                f_error = "output/validation/rec_error_$(col_nm)_g600.nc"
                dist, vals_std = svd_IceSheetDEM.interpolate_raster_to_profile(f_error, xc, yc; band = "std_error")
                vals_std[isnan.(vals_std)] .= 0
                plot!(dist./1e3, vals .- vals_std, fillrange = vals .+ vals_std, label="", fillcolor=color, fillalpha=0.5, lw=0; z_order) #; label, color, ls, lw, ylims, z_order, alpha=0.7)
            end
        end
    end

    # # range of SVD reconstruction from random fields
    # dplot, mid, w = get_range(SVD_RF_files, xc, yc)
    # Plots.plot!(dplot./1e3, mid-w, fillrange=mid+w, fillalpha=0.3, label="SVD RF range", color=cols[3], ls=:dot, linealpha=0)

    # # SEQ simulations
    # dplot, mid, w = get_range(SEQ_files, xc, yc)
    # Plots.plot!(dplot./1e3, mid, label="SGS mean", color=cols[1], ls=:dot, lw=2)
    # Plots.plot!(dplot./1e3, mid-w, fillrange=mid+w, alpha=0.3, label="SGS range", color=cols[1], ls=:dot, linealpha=0)
    push!(ps, p_i)
    Plots.savefig(p_i, joinpath(fig_dir, glacier_name*"_h.png"))
end
Plots.plot(ps..., size=(2300,2000), left_margin = 12Plots.mm, bottom_margin = 12Plots.mm, topmargin = 12Plots.mm)
Plots.savefig(joinpath(fig_dir, "elevation_profiles_all.png"))
# plot two selected glaciers only
i_glaciers = findall(glacier_titles .== "Helheim" .|| glacier_titles .== "Sermeq Kujalleq")
p1_nolegend = plot(ps[i_glaciers[1]], legend=false)
svd_IceSheetDEM.panel_annotate!(p1_nolegend, "b")
p2          = plot(ps[i_glaciers[2]])
svd_IceSheetDEM.panel_annotate!(p2, "a")
plot(p2, p1_nolegend,  wsize=(2000, 800), margin=12Plots.mm)
savefig(joinpath(fig_dir, "elevation_two_glaciers.png"))

# plot dif fields
outline_shp_file = "data/gris-imbie-1980/gris-outline-imbie-1980_updated_crs.shp"
shp              = Shapefile.shapes(Shapefile.Table(outline_shp_file))
f = "output/SVD_reconstruction/dif_lambda_1e7_g600_r700.nc"   # ToDo: make sure this is the latest version with right lambda and r!!
dif = NCDataset(f)["dif"][:,:]
dif[dif .== 0] .= NaN
dif[.!isnan.(dif)] .=  (-1) * dif[.!isnan.(dif)]  # so that it's h_obs - h_SVD
x   = NCDataset(f)["x"][:]
y   = NCDataset(f)["y"][:]

Plots.gr_cbar_width[] = 0.01  # default 0.03

############################################
δx = 80000
δy = 68000

function crds_tup(;yb, xl, δx=80000, δy=68000)
    return (;yb,xl,δx,δy)
end


# coordinates of glacier catchments
coords = Dict("Helheim"         => crds_tup(yb = -2600000, xl =  260000),
              "Sermeq-Kujalleq" => crds_tup(yb = -2310000, xl = -227000),
              "Humboldt"        => crds_tup(yb = -1080000, xl = -380000),
              "79N"             => crds_tup(yb = -1110715, xl =  405378, δy = 110000, δx = 90000),
              "Zachariae"       => crds_tup(yb = -1180715, xl =  435378, δy = 116000, δx = 95000),
              )

function plot_dif(glacier_name, panel_letter, flowline_panel, ins_crds::Tuple, arrow_crds::Tuple)
    iplot        = 1:ixmaxs[glacier_name]
    prof_name    = prof_files[findfirst(occursin.(glacier_name, prof_files))]
    df           = CSV.read(prof_name, DataFrame)

    @unpack yb, xl, δx, δy = coords[glacier_name]
    ix = findall(xl .< x .< xl + δx)
    iy = findall(yb .< y .< yb + δy)
    xtick_interval = 2e4
    xv = Vector(x[ix[1]]:(x[ix[1]]+xtick_interval))
    xtick1 = xv[findfirst(xv .% xtick_interval .== 0)]   # xl+  abs(xtick_interval - (xl .% xtick_interval))
    p_dif = heatmap(x[ix], y[iy], dif[ix,iy]', cmap=:RdBu, clims=(-200,200), aspect_ratio=1, size=(1000,900), margin=10Plots.mm, xlabel="Easting (m)", ylabel="Northing (m)", title=" \n \n"*L"$h_\mathrm{obs}- h_\mathrm{SVD}$", colorbar_title="", grid=false)
    annotate!((xl + 1.25δx, yb + 0.5δy, text("(m)", 18)))
    plot!(p_dif, shp, fill=nothing, xlims=extrema(x[ix]), ylims=extrema(y[iy]), xticks  = xtick1:xtick_interval:x[ix[end]], xtick_direction=:out, lw=0.5)
    Plots.plot!(p_dif, df.X[iplot], df.Y[iplot], label="Flowline in $(flowline_panel))", aspect_ratio=1, lw=3, color=:slategray, legend_foreground_color=nothing, legend_background_color=:transparent) #, ls=:dash)
    svd_IceSheetDEM.panel_annotate!(p_dif, panel_letter)

    rectangle(w, h, x, y) = Shape(x .+ [0,w,w,0], y .+ [0,0,h,h])
    ins = bbox(ins_crds[1], ins_crds[2], 0.2, 0.18, :left)
    plot!(p_dif, inset=ins, subplot=2, aspect_ratio=1)
    plot!(p_dif[2], shp, background_color_inside=nothing, fill=nothing, grid=false, label="", cbar=false, axis=([],false), aspect_ratio=1, lw=0.3)
    plot!(p_dif[2], rectangle(x[ix[end]]-x[ix[1]],y[iy[end]]-y[iy[1]],x[ix[1]],y[iy[1]]), fillalpha=0, linewidth=2, linecolor=:red3, label="")

    # draw arrow manually (didn't manage to adjust arrowhead size with arrow=(..))
    x1, xend = x[ix[1]]+arrow_crds[1], x[ix[1]]+arrow_crds[2]
    y1, yend = y[iy[1]]+arrow_crds[3], y[iy[1]]+arrow_crds[4]
    plot!(p_dif[2], [x1,xend], [y1,yend], color=:red3, lw=2, label="")
    plot!(p_dif[2], [xend,xend], [yend,yend+3e5], color=:red3, lw=2, label="")
    plot!(p_dif[2], [xend,xend+sign(arrow_crds[1]-arrow_crds[2])*3e5], [yend,yend+1e5], color=:red3, lw=2, label="")

    # quiver!(p_dif[2], [x[ix[1]]+arrow_crds[1]],[y[iy[1]]+arrow_crds[3]],quiver=([x[ix[1]]+arrow_crds[2]-x[ix[1]]+arrow_crds[1]],[y[iy[1]]+arrow_crds[4]-y[iy[1]]+arrow_crds[3]]))
    return p_dif
end
p_dif1 = plot_dif("Sermeq-Kujalleq", "c", "a", (0.6,0.6), (5e5, 1.5e5, +6e5, +1.5e5))
# p_dif1 = plot(p_dif1, legend_background_color=:transparent)
# savefig("bla.png")
p_dif2 = plot_dif("Helheim", "d", "b", (0.14,0.6), (-3.5e5, -6e4, +5e5, +1.2e5))
plot(p2, p1_nolegend, p_dif1, p_dif2, wsize=(2500, 1800), left_margin=12Plots.mm) #, top_margin=-15Plots.mm, bottom_margin=10Plots.mm)
savefig(joinpath(fig_dir, "elevation_dif_two_glaciers.png"))


# was a try with bedrock, doesn't show anything particular
# bla = NCDataset("data/bedmachine/BedMachineGreenland-v5.nc")["bed"][:,end:-1:1]
# mask = NCDataset("data/bedmachine/BedMachineGreenland-v5.nc")["mask"][:,end:-1:1]
# xbla = NCDataset("data/bedmachine/BedMachineGreenland-v5.nc")["x"][:]
# ybla = NCDataset("data/bedmachine/BedMachineGreenland-v5.nc")["y"][end:-1:1]
# # ix_bla = findall(coords_Sermeq.xl .< xbla .< coords_Sermeq.xl + δx)
# # iy_bla = findall(coords_Sermeq.yb .< ybla .< coords_Sermeq.yb + δy)
# ix_bla = findall(coords_Helheim.xl .< xbla .< coords_Helheim.xl + δx)
# iy_bla = findall(coords_Helheim.yb .< ybla .< coords_Helheim.yb + δy)
# # prr = heatmap(xbla[ix_bla],ybla[iy_bla], bla[ix_bla,iy_bla]', aspect_ratio=1, cmap=:grays, clims=(-180,500))
# # prr = heatmap(bla[ix_bla,iy_bla]', aspect_ratio=1) #, cmap=:grays, clims=(0,200))
# # plot!(prr, shp, fill=nothing, xlims=extrema(xbla[ix_bla]), ylims=extrema(ybla[iy_bla]), xtick_direction=:out)
# b = nomissing(bla, NaN)
# b[mask .== 0] .= 0.0
# prr = heatmap(xbla[ix_bla],ybla[iy_bla], multihillshade(b[ix_bla,iy_bla], cellsize=150)', aspect_ratio=1, clims=(0, 170), cmap=:grays)
# plot!(prr, shp, fill=nothing, xlims=extrema(xbla[ix_bla]), ylims=extrema(ybla[iy_bla]), xtick_direction=:out)
# plot(prr, p_dif1)



#####################################
# plot map with flowlines           #
#####################################

shp = Shapefile.shapes(Shapefile.Table("data/gris-imbie-1980/gris-outline-imbie-1980_updated_crs.shp"))
p = plot(shp, fillalpha=0, aspect_ratio=1, axis=([],false), wsize = (900, 1400), lw=1) #left_margin = 4Plots.mm, right_margin = 4Plots.mm) #, topmargin = 4Plots.mm)
for fname in prof_files
    df = CSV.read(fname, DataFrame)
    glacier_name = splitext(basename(fname))[1]
    glacier_title = replace(glacier_name, "-" => "\n")
    iplot        = 1:ixmaxs[glacier_name]
    Plots.plot!(df.X[iplot], df.Y[iplot], label="", aspect_ratio=1, lw=3, color="red")
    # Plots.scatter!([df.X[i20[glacier_name]]], [df.Y[i20[glacier_name]]], markersize=4, color="black", label="")
    ann_pos = df.X[iplot[end]] > 0 ? :right : :left
    dif_x   = df.X[iplot[end]] > 0 ? - 5e3 : 2e4
    dif_y   = -2.5e4
    if glacier_name == "Ryder"
        dif_x = 3e4
        dif_y = 0.0
    end
    Plots.annotate!(df.X[iplot[end]]+dif_x, df.Y[iplot[end]]+dif_y, text(glacier_title, ann_pos, 16, "Computer Modern"))
end
Plots.plot(p)
Plots.savefig(joinpath(fig_dir, "map_all.png"))

# smaller insets for examples in main text?? (maybe not)

# p_SK = Plots.plot(p, xlims=(-209467.5,-100027.5), ylims=(-2.31e6,-2.26e6))
# p_HH = Plots.plot(p, xlims=(-209467.5,350027.5), ylims=(-2.65e6,-2.25e6))

# p1 = Plots.heatmap(x[4088:5000],y[9300:9484],vx[4088:5000, 9300:9484]', aspect_ratio=1)
# df = CSV.read("output/profiles/Sermeq-Kujalleq.csv", DataFrame)
# Plots.plot!(df.X,df.Y, label="", lw=2, color="black")
