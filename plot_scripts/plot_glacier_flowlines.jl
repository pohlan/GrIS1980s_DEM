using svd_IceSheetDEM, NCDatasets, Glob, CSV, DataFrames, Shapefile, Plots

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
                "output/SVD_reconstruction/rec_lambda_1e7_g600_r300_dh_detrend.nc",
                "output/geostats_interpolation/kriging/rec_kriging_g600.nc"
               ]
labels = ["Korsgaard et al., 2016", "GrIMP (Howat et al., 2015)", "SVD reconstruction", "Kriging"]
name_for_col =  ["aerodem", "GrIMP", "SVD_dh_detrend", "kriging"]
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
            Plots.scatter!([dist[i20[glacier_name]]./1e3], [vals[i20[glacier_name]]], color="black", markersize=6, label="")
        else
            xmax = Plots.xlims(p_i)[2]
            i_nonan = findall(.!isnan.(vals))
            im   = findmin(abs.(xmax .- dist[i_nonan]./1e3))[2]
            ylims=(-1.0, maximum(vals[i_nonan][1:im])+100)
            Plots.plot!(dist./1e3, vals; label, color, ls, lw, ylims, z_order)
            if label == "Kriging"
                dist, vals_std = svd_IceSheetDEM.interpolate_raster_to_profile(f, xc, yc; band = "std_error")
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
    Plots.savefig(p_i, joinpath(fig_dir, glacier_name*".png"))
end
Plots.plot(ps..., size=(2300,2000), left_margin = 12Plots.mm, bottom_margin = 12Plots.mm, topmargin = 12Plots.mm)
Plots.savefig(joinpath(fig_dir, "elevation_profiles_all.png"))
# plot two selected glaciers only
i_glaciers = findall(glacier_titles .== "Helheim" .|| glacier_titles .== "Sermeq Kujalleq")
p1 = plot(ps[i_glaciers[1]])
svd_IceSheetDEM.panel_annotate!(p1, "a")
p2_nolegend = plot(ps[i_glaciers[2]], legend=false)
svd_IceSheetDEM.panel_annotate!(p2_nolegend, "b")
plot(p1, p2_nolegend,  wsize=(2000, 800), margin=12Plots.mm)
savefig(joinpath(fig_dir, "elevation_two_glaciers.png"))

#####################################
# plot map with flowlines           #
#####################################

shp = Shapefile.shapes(Shapefile.Table("data/gris-imbie-1980/gris-outline-imbie-1980_updated_crs.shp"))
p = heatmap(xx[2500:3900], yy[6400:7900], dh_map[2500:3900, 6400:7900]', clims=(-200,200), aspect_ratio=1, cmap=:coolwarm)
plot!(shp, fillalpha=0, aspect_ratio=1, xlims=(-3e5, -1.5e5), ylims=(-2.3e6, -2.25e6), axis=([],false), wsize = (1600, 400), lw=1) #left_margin = 4Plots.mm, right_margin = 4Plots.mm) #, topmargin = 4Plots.mm)
for fname in prof_files
    df = CSV.read(fname, DataFrame)
    glacier_name = splitext(basename(fname))[1]
    glacier_title = replace(glacier_name, "-" => "\n")
    iplot        = 1:ixmaxs[glacier_name]
    Plots.plot!(df.X[iplot], df.Y[iplot], label="", aspect_ratio=1, lw=3, color="red")
    # Plots.scatter!([df.X[i20[glacier_name]]], [df.Y[i20[glacier_name]]], markersize=4, color="black", label="")
    ann_pos = df.X[iplot[end]] > 0 ? :right : :left
    dif_x   = df.X[iplot[end]] > 0 ? - 5e3 : 2e4
    Plots.annotate!(df.X[iplot[end]]+dif_x, df.Y[iplot[end]]-2.5e4, text(glacier_title, ann_pos, 18, "serif"))
end
Plots.plot(p)
Plots.savefig(joinpath(fig_dir, "map_all.png"))

# smaller insets for examples in main text?? (maybe not)

# p_SK = Plots.plot(p, xlims=(-209467.5,-100027.5), ylims=(-2.31e6,-2.26e6))
# p_HH = Plots.plot(p, xlims=(-209467.5,350027.5), ylims=(-2.65e6,-2.25e6))

# p1 = Plots.heatmap(x[4088:5000],y[9300:9484],vx[4088:5000, 9300:9484]', aspect_ratio=1)
# df = CSV.read("output/profiles/Sermeq-Kujalleq.csv", DataFrame)
# Plots.plot!(df.X,df.Y, label="", lw=2, color="black")
