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
    get_flowline_coords(vx, vy, x0, y0, Δx, Δy, step, L, px0, py0, glacier_name)

### Inputs
- `vx`,`vy`      : matrix of velocities in x,y-direction
- `x0`,`y0`      : first x,y-coordinate of `vx` and `vy` fields
- `Δx`,`Δy`      : cell size in x,y direction of `vx` and `vy`
- `step`         : step length between two points in the flowline
- `L`            : length of profile
- `px0`,`py0`    : coordinates of the starting point for the backpropagating algorithm
- `glacier_name` : name of the glacier

### Output
- String of the csv file the coordinates are saved in
"""
function get_flowline_coords(vx::Matrix, vy::Matrix, x0::T, y0::T, Δx::T, Δy::T, step::T, L::T, px0::T, py0::T, glacier_name::String) where T <: Real
    # set up output directory
    output_dir = joinpath("output", "profiles")
    mkpath(output_dir)
    fname      = joinpath(output_dir, glacier_name*".csv")

    # return file if it exists already
    if isfile(fname)
        return fname
    end

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
        # take a step in opposite direction of velocity vector, with length = step
        push!(px, px[end] - step / v_norm * vx_interp)
        push!(py, py[end] - step / v_norm * vy_interp)
        nit += 1
        dist_tot += step
    end
    # save as csv
    df_coords = DataFrame("X" => px, "Y" => py)
    CSV.write(fname, df_coords)
    return fname
end

"""
    extract_profile_nearest_points(fname, xc, yc; band=1) -> (dist,vals)

Find points on raster corresponding to flowline profile; no interpolation, just assigns values of nearest neighbor cells.
### Inputs
- `fname`   : filename of raster field that will be profiled
- `xc`,`yc` : coordinates of the flowline profile
- `band`    : band of the field in `fname`, usually 1 if there's only one band

### Outputs
- `dist`    : vector of distances along profile
- `vals`    : field values corresponding to these points along the profile
"""
function extract_profile_nearest_points(fname, xc, yc; band=1)   # inspired by https://github.com/evetion/GeoArrays.jl/blob/master/src/geoutils.jl#L270
    # read in target field as GeoArray
    ga     = GeoArrays.read(fname, masked=false)

    # for each point of flowline, find indices of closest raster cell
    coor   = [GeoArrays.indices(ga, (xc[i],yc[i])).I for i in eachindex(xc)]
    unique!(coor)

    # loop through pairs of subsequent points along flowline
    vals = Float32[]
    dist   = Float32[]
    dist_0 = 0
    for (start_pt, end_pt) in zip(coor[1:end-1], coor[2:end])
        i0, j0 = start_pt      # indices in x-direction of two subsequent raster cells along profile
        i1, j1 = end_pt        # indices in y-direction ""
        δx     = i1 - i0
        δy     = j1 - j0
        pos0   = GeoArrays.coords(ga, (i0, j0))
        if abs(δx) >= abs(δy)
            j = j0
            xstep = δx > 0 ? 1 : -1
            for (d, i) in enumerate(i0:xstep:i1)
                idx = i, j+div((d - 1) * δy, abs(δx), RoundNearest), band
                pos = GeoArrays.coords(ga, (idx[1:2]))
                tot_dist = norm(pos .- pos0)+dist_0
                push!(vals, ga[idx...])
                push!(dist, tot_dist)
            end
        else
            i = i0
            ystep = δy > 0 ? 1 : -1
            for (d, j) in enumerate(j0:ystep:j1)
                idx = i+div((d - 1) * δx, abs(δy), RoundNearest), j, band
                pos = GeoArrays.coords(ga, (idx[1:2]))
                tot_dist = norm(pos .- pos0)+dist_0
                push!(vals, ga[idx...])
                push!(dist, tot_dist)
            end
        end
        dist_0 = dist[end]
    end
    return dist, vals
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
    Z = replace_missing(Z, NaN)
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
