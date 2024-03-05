using svd_IceSheetDEM
using NCDatasets, Interpolations, DataFrames, CSV, ProgressMeter, GeoStats, LsqFit, ParallelStencil
using StatsBase
import Plots

@init_parallel_stencil(Threads, Float64, 1)

gr = 1200

get_ix(i,nx) = i % nx == 0 ? nx : i % nx
get_iy(i,nx) = cld(i,nx)


# masks
ds_mask = NCDataset("data/gris-imbie-1980/imbie_mask_g$(gr).nc")["Band1"][:]
bedm_mask = NCDataset("data/bedmachine/bedmachine_g$(gr).nc")["mask"][:]

# GrIMP
grimp = NCDataset("data/bedmachine/bedmachine_g$(gr).nc")["surface"][:,:]

# aero
ds_aero = NCDataset("data/grimp_minus_aero_g$(gr).nc")
dh_aero_all = ds_aero["dh"][:,:,1]
dh_aero_all = dh_aero_all[:,end:-1:1]
x  = sort(ds_aero["x"])
y  = sort(ds_aero["y"])
idx_aero = findall(.!ismissing.(vec(dh_aero_all)) .&& .!ismissing.(vec(ds_mask)) .&& vec((bedm_mask .!= 1)))
df_aero = DataFrame(:x       => x[get_ix.(idx_aero, length(x))],
                        :y       => y[get_iy.(idx_aero, length(x))],
                        :h_grimp => grimp[idx_aero],
                        :dh      => dh_aero_all[idx_aero],
                        :idx     => idx_aero)
# ii = sort(unique(rand(1:length(idx_aero), fld(size(df_aero,1), 2))))
# keepat!(df_aero, ii)
# remove outliers              # if outliers are removed the interpolated points will have significant differences to original aerodem data
# aero_to_delete = findall(abs.(df_aero.dh) .> 7 .* mad(df_aero.dh)) # .|| df_aero.h_grimp .< 50)
# deleteat!(df_aero, aero_to_delete)

# atm
df_atm = CSV.read("atm_grimp_interp_pts.csv", DataFrame)
mindist = 5e4 # minimum distance between aerodem and atm points
function choose_atm(x, y)::Bool
    dist_to_aero = maximum(pairwise(Distances.euclidean, [x y], [df_aero.x df_aero.y], dims=1)[:])
    is_above_mindist = dist_to_aero > mindist
    return is_above_mindist
end
println("selecting flightline values...")
id = unique(sort(rand(1:size(df_atm,1), 100000)))
keepat!(df_atm, id)
filter!([:x,:y] => choose_atm, df_atm)
# remove outliers
atm_to_delete = findall(abs.(df_atm.dh) .> 5 .* mad(df_atm.dh))
deleteat!(df_atm, atm_to_delete)

# merge aerodem and atm data
df_all = vcat(df_aero, df_atm, cols=:intersect)

# plot
dh_plot = Array{Union{Float32,Missing}}(missing, (length(x),length(y)))
dh_plot[df_aero.idx] = df_aero.dh
Plots.heatmap(x,y, dh_plot', cmap=:bwr, clims=(-10,10))
Plots.scatter!(df_atm.x, df_atm.y, marker_z=df_atm.dh, color=:bwr, markersize=0.5, markerstrokewidth=0, legend=false)
Plots.savefig("dplot.png")

# mean trend as a fct of elevation
dh_binned, ELV_bin_centers = svd_IceSheetDEM.bin_equal_sample_size(df_all.h_grimp, df_all.dh, 12000)  # 16000
itp_bias = linear_interpolation(ELV_bin_centers, mean.(dh_binned),  extrapolation_bc = Interpolations.Flat())
Plots.scatter(ELV_bin_centers, mean.(dh_binned))
ELV_plot = minimum(df_all.h_grimp):10:maximum(df_all.h_grimp)
Plots.plot!(ELV_plot, itp_bias(ELV_plot), label="linear interpolation", legend=:bottomright)
Plots.savefig("c1plot.png")

# standard deviation as a fct of elevation
itp_std = linear_interpolation(ELV_bin_centers, std.(dh_binned),  extrapolation_bc = Interpolations.Flat())
Plots.scatter(ELV_bin_centers, std.(dh_binned))
ELV_plot = minimum(df_all.h_grimp):10:maximum(df_all.h_grimp)
Plots.plot!(ELV_plot, itp_std(ELV_plot), label="linear interpolation", legend=:topright)
Plots.savefig("c2plot.png")

# standardize
df_all.dh_detrend = (df_all.dh .- itp_bias(df_all.h_grimp)) ./ itp_std(df_all.h_grimp)
Plots.density(df_all.dh_detrend, label="Standardized observations")
Plots.density!((df_all.dh .- mean(df_all.dh)) ./ std(df_all.dh), label="Standardized without binning")
Plots.plot!(Normal(), label="Normal distribution")
Plots.savefig("eplot.png")

# variogram
println("Variogram...")
table_all  = (; Z=df_all.dh_detrend)
coords_all = [(df_all.x[i], df_all.y[i]) for i in 1:length(df_all.x)]
data       = georef(table_all,coords_all)
gamma      = EmpiricalVariogram(data, :Z; nlags=400,  maxlag=2e5)
# fit a covariance function
# function custom_var(x, params)
#     γ1 = ExponentialVariogram(range=params[1], sill=params[2], nugget=params[3])      # forcing nugget to be zero, not the case with the built-in fit function of GeoStats
#     # γ2 = ExponentialVariogram(range=params[3], sill=params[4], nugget=0.0)
#     f = γ1.(x) #.+ γ2.(x)
#     return f
# end
# ff = LsqFit.curve_fit(custom_var, gamma.abscissa, gamma.ordinate, [5e4, 1.0, 0.2]);
# julia> ff.param
# 2-element Vector{Float64}:
#  50280.021902654196
#      0.94295155705642
varg = ExponentialVariogram(range=ff.param[1], sill=ff.param[2], nugget=ff.param[3])
varg_gs_exp = fit(ExponentialVariogram, gamma)
varg_gs_any = fit(Variogram, gamma)
Plots.scatter(gamma.abscissa, gamma.ordinate, label="observations", xscale=:log10)
# Plots.plot!(gamma.abscissa,varg.(gamma.abscissa), label="LsqFit fit")
Plots.plot!(gamma.abscissa,varg_gs_exp.(gamma.abscissa), label="GeoStats fit exponential")
Plots.plot!(gamma.abscissa,varg_gs_any.(gamma.abscissa), label="GeoStats fit circular")
Plots.savefig("vplot.png")

# varg = ExponentialVariogram(range=50280.021902654196, sill=0.94295155705642, nugget=0.05)   # nugget is uncertainty of observations
varg = varg_gs_any

@parallel_indices (ibox) function interpolate_subsets!(mb, dh_boxes, x_boxes, y_boxes, x_points, y_points, gr)
    if 1 <= ibox <= length(x_boxes)
        # extract data points in sub-area
        x_min, x_max = extrema(x_points[ibox])
        y_min, y_max = extrema(y_points[ibox])

        # prepare data and grid for interpolation
        table_box  = (; Z=Float32.(dh_boxes[ibox]))
        coords_box = [(xi,yi) for (xi,yi) in zip(x_boxes[ibox],  y_boxes[ibox])]
        geotable   = georef(table_box,coords_box)
        grid       = CartesianGrid((x_min-gr/2, y_min-gr/2), (x_max+gr/2, y_max+gr/2), (Float64(gr), Float64(gr)))
        model      = Kriging(varg_gs_any, 0.0)

        # do interpolation
        interp = geotable |> Interpolate(grid, model) #, minneighbors=300, maxneighbors=500)
        mb[ibox,:,:] = reshape(interp.Z, size(grid))
    end
    return
end

# x_i1 = 110:450
# y_i1 = 1550:1800
# # 3.6 hours

# mb_1 = interpolate_subset(df_all, x_i1, y_i1)

# x_i2 = 170:510
# y_i2 = 1550:1800

# mb_2 = interpolate_subset(df_all, x_i2, y_i2)

x_is = [110:280, 170:340, 230:400, 290:460, 350:520, 410:580, 470:640, 530:700]
y_is = [1550:1670, 1550:1670, 1550:1670, 1550:1670, 1550:1670, 1550:1670, 1550:1670, 1550:1670]

x_boxes = []
y_boxes = []
x_points = []
y_points = []
dh_boxes = []
for (ix, iy) in zip(x_is, y_is)
    x_min, x_max = extrema(x[ix])
    y_min, y_max = extrema(y[iy])
    function is_in_box(x, y)::Bool
        is_in_x = x_min < x < x_max
        is_in_y = y_min < y < y_max
        return is_in_x && is_in_y
    end
    df_box = filter([:x,:y] => is_in_box, df_all)
    # df_aero_box = filter([:x,:y] => is_in_box, df_aero)
    # df_aero_box.dh_detrend = (df_aero_box.dh .- itp_bias(df_aero_box.h_grimp)) ./ itp_std(df_aero_box.h_grimp)
    # iix = [findfirst(i .== x[x_i]) for i in df_aero_box.x]
    # iiy = [findfirst(i .== y[y_i]) for i in df_aero_box.y]
    push!(x_boxes, df_box.x)
    push!(y_boxes, df_box.y)
    push!(dh_boxes, df_box.dh_detrend)
    push!(x_points, x[ix])
    push!(y_points, y[iy])
end



mb   = @zeros(length(x_boxes), length(x_is[1]), length(y_is[1]))
tic = Base.time()
@parallel interpolate_subsets!(mb, dh_boxes, x_boxes, y_boxes, x_points, y_points, gr)
toc = Base.time() - tic
tt = toc / 60
println("Interpolation took $tt minutes.")





## Plotting stuff


# Plots.heatmap(x[x_i], y[y_i], mb', cmap=:bwr, clims=(-5,5), aspect_ratio=1)


overlap_1 = mb_1[61:end,:]
overlap_2 = mb_2[1:end-60,:]
Plots.heatmap(x[x_i1[61:end]], y[y_i1],overlap_1' .- overlap_2', cmap=:bwr, clims=(-0.5,0.5))
Plots.savefig("diff_overlap.png")

# geotable, grid, model, xis, yis, grimp, itp_bias, itp_std, idx_aero = prepare(gr)
# nx, ny = length(xis), length(yis)
# mb = do_kriging(geotable, grid, model, nx,ny)

h_predict = grimp[x_i,y_i] .- (mb.*itp_std.(grimp[x_i,y_i]) .+ itp_bias.(grimp[x_i,y_i]))
Plots.heatmap(x[x_i], y[y_i], h_predict', aspect_ratio=1, clims=(0,1600))
# Plots.savefig("h_predict.png")

h_aero = grimp .- dh_aero_all
idx_box = findall(.!ismissing.(h_aero[x_i,y_i]))


dd = Array{Union{Float32,Missing}}(missing, (length(x_i),length(y_i)))
dd[idx_box_dhdata] = h_predict[idx_box_dhdata] .- h_aero[x_i,y_i][idx_box_dhdata]
Plots.heatmap(x[x_i], y[y_i], dd', aspect_ratio=1, cmap=:bwr, clims=(-60,60))
Plots.savefig("diff_to_aero.png")

Plots.heatmap(x[x_i],y[y_i],h_aero[x_i,y_i]', aspect_ratio=1, clims=(0,1500))
hh = Array{Union{Float32,Missing}}(missing, (length(x_i),length(y_i)))
hh[idx_box] = h_predict[idx_box]
Plots.heatmap(x[x_i],y[y_i],hh', aspect_ratio=1, clims=(0,1500))


# x_min, x_max = extrema(x[x_i])
# y_min, y_max = extrema(y[y_i])
# df_aero_box = filter([:x,:y] => is_in_box, df_aero)
db = Array{Union{Float32,Missing}}(missing, (length(x),length(y)))
db[df_aero_box.idx] .= df_aero_box.dh
Plots.heatmap(x[x_i], y[y_i],reshape(db,length(x),length(y))[x_i,y_i]', cmap=:bwr, clims=(-5,5), aspect_ratio=1)
idx_box_dhdata = findall(.!ismissing.(reshape(db,length(x),length(y))[x_i,y_i]))

mbpart = zeros(length(x_i),length(y_i))
mbpart[idx_box_dhdata] .= mb[idx_box_dhdata]
Plots.heatmap(x[x_i], y[y_i],mbpart', cmap=:bwr, clims=(-5,5), aspect_ratio=1)
