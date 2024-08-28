using svd_IceSheetDEM, NCDatasets, Glob, CSV, DataFrames, Shapefile

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
                       "Sermeq-Kujalleq"         => (-197941., -2271719.),
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
                "output/rec_files/rec_lambda_1e7_g600_r200.nc",
                "output/geostats_interpolation/kriging/rec_kriging_g600.nc"
               ]
labels = ["Korsgaard et al., 2016", "GrIMP (Howat et al., 2015)", "SVD reconstruction", "Kriging"]
bandnm = ["Band1", "surface", "surface", "surface"]
cols   = [:darkorchid, :orange, :olive, :orchid]
lstls  = [:solid, :solid, :dot, :dot]
lw     = 4
z_orders = [1,1,1,1]

Plots.scalefontsizes()
Plots.scalefontsizes(1.9)
ps = Plots.Plot{Plots.GRBackend}[]

ixmaxs = Dict()
for (ip, pf) in enumerate(prof_files) #[[3,4,5,6]])
    glacier_name = splitext(basename(pf))[1]
    glacier_title = replace(glacier_name, "-" => " ")
    prof = CSV.read(pf, DataFrame)
    xc = prof.X
    yc = prof.Y

    # set up figure
    xlabel = "Distance along profile [km]"
    ylabel = "Surface elevation [m]"
    legend = :topleft
    p_i = Plots.plot(title="\n"*glacier_title, size=(1200,1000); xlabel, ylabel, left_margin = 12Plots.mm, bottom_margin = 4Plots.mm, topmargin = 4Plots.mm,
              fg_legend=:transparent, bg_legend=:transparent, legend)

    # plot GrIMP, reconstruction and aerodem
    for (i,(f, label, color, ls, z_order, band)) in enumerate(zip(files, labels, cols, lstls, z_orders, bandnm))
        dist, vals   = svd_IceSheetDEM.interpolate_raster_to_profile(f, xc, yc; band)
        if label == "Korsgaard et al., 2016"
            i_nonan = findall(.!isnan.(vals))
            i_max   = findmax(vals[i_nonan])[2]
            xmax    = min(dist[i_nonan[i_max]]* 1.5, maximum(dist))
            ixmaxs[glacier_name] = i_nonan[i_max]
            xlims=(0, xmax/1e3)
            Plots.plot!(dist./1e3, vals; label, color, ls, lw, z_order, xlims)
        else
            xmax = Plots.xlims(p_i)[2]
            im   = findmin(abs.(xmax .- dist./1e3))[2]
            ylims=(-1.0, maximum(vals[1:im])+200)
            Plots.plot!(dist./1e3, vals; label, color, ls, lw, ylims, z_order)
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
Plots.plot(ps..., size=(2000,1900), left_margin = 12Plots.mm, bottom_margin = 12Plots.mm, topmargin = 12Plots.mm)
Plots.savefig(joinpath(fig_dir, "elevation_profiles_all.png"))

#####################################
# plot map with flowlines           #
#####################################

shp = Shapefile.shapes(Shapefile.Table("data/gris-imbie-1980/gris-outline-imbie-1980_updated_crs.shp"))
p = plot(shp, fillalpha=0, aspect_ratio=1, xlims=(-7e5, 9e5), axis=([],false), wsize = (800, 1300), lw=1) #left_margin = 4Plots.mm, right_margin = 4Plots.mm) #, topmargin = 4Plots.mm)
for fname in prof_files
    df = CSV.read(fname, DataFrame)
    glacier_name = splitext(basename(fname))[1]
    glacier_title = replace(glacier_name, "-" => "\n")
    iplot        = 1:ixmaxs[glacier_name]
    Plots.plot!(df.X[iplot], df.Y[iplot], label="", aspect_ratio=1, lw=3, color="red")
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
