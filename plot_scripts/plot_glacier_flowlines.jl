using svd_IceSheetDEM, NCDatasets, Glob, CSV, DataFrames

#####################################
# calculate flowlines from velocity #
#####################################

# read in
vx = NCDataset("data/ITS_LIVE_velocity_120m_RGI05A_0000_v02_vx.nc")["vx"][:,:]
vy = NCDataset("data/ITS_LIVE_velocity_120m_RGI05A_0000_v02_vy.nc")["vy"][:,:]

x  = NCDataset("data/ITS_LIVE_velocity_120m_RGI05A_0000_v02_vx.nc")["x"][:]
x0 = x[1]
Δx = x[2]-x[1]
y  = NCDataset("data/ITS_LIVE_velocity_120m_RGI05A_0000_v02_vx.nc")["y"][:]
y0 = y[1]
Δy = y[2]-y[1]


step = 200.0  # m
px0, py0 = x[4130], y[9417]
maxdist = 8e4

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
    fname = svd_IceSheetDEM.get_flowline_coords(vx, vy, x0, y0, Δx, Δy, step, maxdist, px0, py0, output_name)
    push!(prof_files, fname)
end

# plot
p = Plots.plot()
for fname in prof_files
    df = CSV.read(fname, DataFrame)
    print(fname*": "); println(length(df.X))
    Plots.plot!(df.X, df.Y, label="", aspect_ratio=1, lw=2, color="red")
end
Plots.plot(p)

p1 = Plots.heatmap(x[4088:4567],y[9300:9484],vx[4088:4567, 9300:9484]')
df = CSV.read("output/profiles/Sermeq-Kujalleq.csv", DataFrame)
Plots.plot!(df.X,df.Y, label="", lw=2, color="white")

#####################################
# plot elevations along flowlines   #
#####################################


# randn_files = glob("output/rec_rand_id_*.nc")
SEQ_files    = glob("output/SEQ/simulations/SEQ_maxngh_100_g600_id_*_upsampled.nc")
SVD_RF_files = glob("output/SVD_RF/simulations/rec_rand_id_*.nc")
files        = ["data/aerodem/aerodem_rm-filtered_geoid-corr_g600.nc",
                "data/bedmachine/bedmachine_g600.nc",
                "output/rec_files/rec_lambda_1e7_g600_r200.nc",
                "output/geostats_interpolation/kriging/rec_kriging_g600.nc"
               ]
labels = ["Korsgaard et al., 2016", "GrIMP (Howat et al., 2015)", "SVD reconstruction", "Kriging"]
bandnm = ["Band1", "surface", "surface", "surface"]
cols   = [:black, :orange, :teal, :orchid]
lstls  = [:solid, :solid, :dot, :dot]
lw     = 3

Plots.scalefontsizes()
Plots.scalefontsizes(1.9)
ps = Plots.Plot{Plots.GRBackend}[]
for (ip, pf) in enumerate(prof_files) #[[3,4,5,6]])
    glacier_name = split(basename(pf), "_")[1]
    prof = CSV.read(pf, DataFrame)
    xc = prof.X
    yc = prof.Y

    # set up figure
    xlabel = (ip == 3 || ip == 4) ? "Distance along profile [km]" : ""
    ylabel = isodd(ip)  ? "Surface elevation [m]" : ""
    legend = ip == 1    ? :topleft : false
    p_i = Plots.plot(title="\n"*splitext(glacier_name)[1], size=(1000,1000); xlabel, ylabel, left_margin = 12Plots.mm, bottom_margin = 4Plots.mm, topmargin = 4Plots.mm,
              fg_legend=:transparent, bg_legend=:transparent, legend)

    # plot GrIMP, reconstruction and aerodem
    for (i,(f, label, color, ls, band)) in enumerate(zip(files, labels, cols, lstls, bandnm))
        dist, vals   = svd_IceSheetDEM.interpolate_raster_to_profile(f, xc, yc; band)
        # i_to_plot    = findall(vals.>0)
        if label == "Korsgaard et al., 2016"
            # xlims=(0, min(maximum(dist[i_to_plot])/1e3 * 1.5, 45))
            Plots.plot!(dist./1e3, vals; label, color, ls, lw) #, xlims)
        else
            xmax = Plots.xlims(p_i)[2]
            # im   = findmin(abs.(xmax .- dist./1e3))[2]
            # ylims=(0, vals[i_to_plot[im]])
            Plots.plot!(dist./1e3, vals; label, color, ls, lw) #, ylims)
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
end

Plots.plot(ps..., size=(1600,1100), left_margin = 12Plots.mm, bottom_margin = 12Plots.mm, topmargin = 12Plots.mm)
Plots.savefig("flowline_profiles.png")
