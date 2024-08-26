
"""
    get_flowline_coords(vx, vy, x0, y0, Δx, Δy, step, L, px0, py0, glacier_name)

### Inputs:
- `vx`,`vy`      : matrix of velocities in x,y-direction
- `x0`,`y0`      : first x,y-coordinate of `vx` and `vy` fields
- `Δx`,`Δy`      : cell size in x,y direction of `vx` and `vy`
- `step`         : step length between two points in the flowline
- `L`            : length of profile
- `px0`,`py0`    : coordinates of the starting point for the backpropagating algorithm
- `glacier_name` : name of the glacier

### Output:
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

    # define a function to interpolate the velocity to any point, depends on x0, y0 and gr
    function interp_iwd(v::Matrix, px, py)
        # calculate horizontal/vertical distance to left x points and lower y points
        xl = mod(px-x0, Δx)   # horizontal distance to left points
        yl = mod(py-y0, Δy)   # vertical distance to lower points

        # calculate weights for all neighboring points
        wt_ll = (1-xl/Δx) * (1-yl/Δy)  # lower left
        wt_ur =    xl/Δx  *    yl/Δy   # upper right
        wt_ul = (1-xl/Δx) *    yl/Δy   # upper left
        wt_lr =    xl/Δx  * (1-yl/Δy)  # lower right

        # determine corner indices and velocities
        ixl = Int(fld(px-x0, Δx)) + 1   # x indices of left points
        iyl = Int(fld(py-y0, Δy)) + 1   # y indices of right points

        # calculate velocity based on weights
        sum_v_wt = v[ixl,iyl]*wt_ll + v[ixl+1,iyl+1]*wt_ur + v[ixl,iyl+1]*wt_ul + v[ixl+1,iyl]*wt_lr
        sum_wt   = wt_ll + wt_ur + wt_ul + wt_lr
        v_interp = sum_v_wt / sum_wt
        return v_interp
    end

    # backtracking of velocity to derive flow profile
    px, py   = [px0], [py0]
    nit      = 0
    maxit    = 500
    dist_tot = 0
    while dist_tot < L && nit < maxit
        vx_interp = interp_iwd(vx, px[end], py[end])
        vy_interp = interp_iwd(vy, px[end], py[end])
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
