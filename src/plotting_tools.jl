# Define colors for different DEMs / methods for DEM reconstruction
pt1 = Plots.palette(:batlow10)
pt2 = Plots.palette(:bamO25)
# pt2 = Plots.palette(:bamO25)
pt3 = Plots.palette(:berlin10)
palette_dict = Dict("combined"       => pt2[6],
                    "GP"             => pt3[2],
                    "SVD_dh_detrend" => pt2[2],
                    "SVD_dh"         => pt2[4],
                    "SVD"            => pt2[6],
                    "aerodem"        => pt1[1],
                    "GrIMP"          => pt1[7])

font_scaling = 1.9
wwidth  = 1000
wheight = 700

panel_annotate!(p, letter) = annotate!(p, (xlims(p)[1]+(xlims(p)[2]-xlims(p)[1])*0.02, ylims(p)[1]+(ylims(p)[2]-ylims(p)[1])*1.07, Plots.text(L"\textbf{(%$letter)}", :left, p.attr[:plot_titlefontsize]-4)))
panel_annotate_xlog!(p, letter) = annotate!(p, (10 .^(log10(xlims(p)[1])+(log10(xlims(p)[2])-log10(xlims(p)[1]))*0.02), ylims(p)[1]+(ylims(p)[2]-ylims(p)[1])*1.07, Plots.text(L"\textbf{(%$letter)}", :left, p.attr[:plot_titlefontsize]-4)))
panel_annotate_ylog!(p, letter) = annotate!(p, (xlims(p)[1]+(xlims(p)[2]-xlims(p)[1])*0.02, 10 .^(log10(ylims(p)[1])+(log10(ylims(p)[2])-log10(ylims(p)[1]))*1.07), Plots.text(L"\textbf{(%$letter)}", :left, p.attr[:plot_titlefontsize]-4)))

"""
    interp_raster(A, px, py, x0, y0, Δx, Δy) -> A_interp

Interpolation of raster `A` onto an arbitrary point (`px`,`py`) in the domain.
Finds the four neighboring raster points and calculates a weighted average:

A_interp = ∑(A[i,j]*wt[i,j]) / ∑(wt[i,j])

where wt[i,j] are weights based on distance to each point.

### Inputs
- `A`       : raster given as a matrix
- `px`,`py` : coordinates of point to interpolate on
- `x0`,`y0` : first x,y coordinates of raster
- `Δx`,`Δy` : cell size in x,y direction of the raster

### Output
- `A_interp` : interpolated value of `A` in point (`px`,`py`)
"""
function interp_raster(A::Matrix, px::T, py::T, x0::T, y0::T, Δx::T, Δy::T) where T <: Real
    # calculate horizontal/vertical distance to left x points and lower y points
    xl = mod(px-x0, Δx)   # horizontal distance to left points
    yl = mod(py-y0, Δy)   # vertical distance to lower points

    # calculate weights for all neighboring points
    wt_ll = (1-xl/Δx) * (1-yl/Δy)  # lower left
    wt_ur =    xl/Δx  *    yl/Δy   # upper right
    wt_ul = (1-xl/Δx) *    yl/Δy   # upper left
    wt_lr =    xl/Δx  * (1-yl/Δy)  # lower right

    # determine corner indices
    ixl = Int(fld(px-x0, Δx)) + 1   # x indices of left points
    iyl = Int(fld(py-y0, Δy)) + 1   # y indices of right points

    # calculate field values based on weights
    sum_A_wt = A[ixl,iyl]*wt_ll + A[ixl+1,iyl+1]*wt_ur + A[ixl,iyl+1]*wt_ul + A[ixl+1,iyl]*wt_lr
    sum_wt   = wt_ll + wt_ur + wt_ul + wt_lr
    A_interp = sum_A_wt / sum_wt
    return A_interp
end

"""
    get_flowline_coords(vx, vy, x0, y0, Δx, Δy, lstep, L, px0, py0, glacier_name)

### Inputs
- `vx`,`vy`      : matrix of velocities in x,y-direction
- `x0`,`y0`      : first x,y-coordinate of `vx` and `vy` fields
- `Δx`,`Δy`      : cell size in x,y direction of `vx` and `vy`
- `lstep`         : lstep length between two points in the flowline
- `L`            : length of profile
- `px0`,`py0`    : coordinates of the starting point for the backpropagating algorithm
- `glacier_name` : name of the glacier

### Output
- String of the csv file the coordinates are saved in
"""
function get_flowline_coords(vx_file::String, vy_file::String, lstep::T, L::T, px0::T, py0::T, glacier_name::String) where T <: Real
    # set up output directory
    output_dir = joinpath("output", "profiles")
    mkpath(output_dir)
    fname      = joinpath(output_dir, glacier_name*".csv")

    # return file if it exists already
    if isfile(fname)
        return fname
    end

    # read in
    vx = NCDataset(vx_file)["vx"][:,:]
    vy = NCDataset(vy_file)["vy"][:,:]
    x  = NCDataset(vx_file)["x"][:]
    x0 = x[1]
    Δx = x[2]-x[1]
    y  = NCDataset(vx_file)["y"][:]
    y0 = y[1]
    Δy = y[2]-y[1]

    # define a function to interpolate the velocity to any point, depends on x0, y0, Δx and Δy
    interp_v(v, px, py) = interp_raster(v, px, py, x0, y0, Δx, Δy)

    # backtracking of velocity to derive flow profile
    px, py   = [px0], [py0]
    nit      = 0
    maxit    = 500
    dist_tot = 0
    while dist_tot < L && nit < maxit
        vx_interp = interp_v(vx, px[end], py[end])
        vy_interp = interp_v(vy, px[end], py[end])
        v_norm    = sqrt(vx_interp^2 + vy_interp^2)
        if any(ismissing.([vx_interp,vy_interp]))
            break
        end
        # take a step in opposite direction of velocity vector, with length = lstep
        push!(px, px[end] - lstep / v_norm * vx_interp)
        push!(py, py[end] - lstep / v_norm * vy_interp)
        nit += 1
        dist_tot += lstep
    end
    # save as csv
    df_coords = DataFrame("X" => px, "Y" => py)
    CSV.write(fname, df_coords)
    return fname
end

"""
    interpolate_raster_to_profile(fname, xc, yc; band="Band1")

Projects values of raster field onto coordinates of arbitrary profile; for interpolation see `interp_raster()`

### Inputs
- `fname`   : filename (netcdf) of raster
- `xc`,`yc` : coordinates of flowline
- `band`    : name of band in netcdf file

### Outputs
- `dists`     : vector of distances along profile
- `vals`    : field values corresponding to these points along the profile
"""
function interpolate_raster_to_profile(fname, xc, yc; band="Band1")
    # read in netcdf file
    ds = NCDataset(fname)
    Z  = ds[band][:,:]
    Z  = nomissing(Z, NaN)
    Z[Z .== 0.0] .= NaN
    x  = ds["x"][:]
    y  = ds["y"][:]
    x0, y0 = x[1], y[2]
    Δx, Δy = x[2]-x[1], y[2]-y[1]

    interp_z(px, py) = interp_raster(Z, px, py, x0, y0, Δx, Δy)
    dists = zeros(length(xc))
    vals  = zeros(length(xc))
    for (ip, (px,py)) in enumerate(zip(xc,yc))
        vals[ip]  = interp_z(px, py)
        dists[ip] = norm([px-xc[1], py-yc[1]])
    end
    return dists, vals
end
