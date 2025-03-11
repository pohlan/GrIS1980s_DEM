function make_geotable(input::AbstractVector, x::AbstractVector, y::AbstractVector)
    table    = (;Z = input)
    crds   = [(xi, yi) for (xi,yi) in zip(x,y)]
    return georef(table, crds)
end

"""
Binning 1D based on equal bin size
"""
function bin_equal_bin_size(x, y, n_bins)   # 1D
    bin_edges = range(extrema(x)..., n_bins+1)
    out = Vector{Vector{eltype(y)}}(undef,n_bins)
    for b in 1:n_bins
        i = findall(bin_edges[b] .< x .< bin_edges[b+1])
        out[b] = y[i]
    end
    bin_centers = collect(0.5.*(bin_edges[1:end-1] .+ bin_edges[2:end]))

    id_del = findall(length.(out) .< 100)
    deleteat!(out, id_del)
    deleteat!(bin_centers, id_del)
    return out, bin_centers
end

"""
Binning 1D or 2D based on equal sample size per bin
"""
function bin_equal_sample_size(x, y, n_samples)   # 1D
    nx = length(x)
    n_bins = ceil(Int, nx / n_samples)
    if nx % n_samples < 100
        n_bins -= 1
    end
    p = sortperm(x)
    out = Vector{Vector{eltype(y)}}(undef,n_bins)
    bin_centers = zeros(eltype(x), n_bins)
    for b in 1:n_bins
        if b == n_bins
            i = (b-1)*n_samples+1:nx
        else
            i = (b-1)*n_samples+1:b*n_samples
        end
        bin_centers[b] = 0.5*(sort(x)[i[1]]+sort(x)[i[end]])
        out[b] = y[p[i]]
    end
    return out, bin_centers
end
function bin_equal_sample_size(x1, x2, y, n_bins_1, n_bins_2) # 2D
    @assert length(x1) == length(x2) == length(y)
    nx = length(x1)
    nsampl1 = ceil(Int, nx / n_bins_1)
    nsampl2 = ceil(Int, nx / n_bins_2)
    if nx % nsampl1 < 50
        n_bins_1 -= 1
    end
    if nx % nsampl2 < 50
        n_bins_2 -= 1
    end
    p1 = sortperm(x1)
    p2 = sortperm(x2)
    bin_edges_1 = x1[p1[[1:nsampl1:(n_bins_1-1)*nsampl1+1;nx]]]
    bin_edges_2 = x2[p2[[1:nsampl2:(n_bins_2-1)*nsampl2+1;nx]]]
    out = Array{Vector{eltype(y)}}(undef,n_bins_1,n_bins_2)
    for b1 in 1:n_bins_1
        for b2 in 1:n_bins_2
            i = findall(bin_edges_1[b1] .< x1 .< bin_edges_1[b1+1] .&& bin_edges_2[b2] .< x2 .< bin_edges_2[b2+1])
            out[b1,b2] = y[i]
        end
    end
    bin_centers_1 = bin_edges_1[1:end-1] .+ 0.5*diff(bin_edges_1)
    bin_centers_2 = bin_edges_2[1:end-1] .+ 0.5*diff(bin_edges_2)
    return out, bin_centers_1, bin_centers_2
end

function get_ATM_df(fname, x, y, bin_field_1, df_aero; mindist=1e4, I_no_ocean)
    df_atm = CSV.read(fname, DataFrame)
    df_atm[!,:source] .= :atm
    df_atm[!,:h]      .= df_atm.h_ref .- df_atm.dh
    sort!(unique!(x))   # necessary for interpolation
    sort!(unique!(y))
    itp = interpolate((x, y), bin_field_1, Gridded(Linear()))
    df_atm[!,:bfield_1] = itp.(df_atm.x, df_atm.y)
    df_atm[!,:bfield_2] = df_atm.h_ref
    # prepare neighbor search to filter out values too close to aerodem or overlapping
    d        = PointSet([Point(xi, yi) for (xi,yi) in zip(df_aero.x,df_aero.y)])
    searcher = BallSearch(d, MetricBall(mindist))
    function choose_atm(x_, y_)::Bool
        ix    = findmin(abs.(x_ .- x))[2]
        iy    = findmin(abs.(y_ .- y))[2]
        iglob = get_global_i(ix, iy, length(x))
        # filter values outside the ice sheet
        is_in_icesheet = iglob ∈ I_no_ocean
        # filter values overlapping or close to aerodem
        inds = search(Point(x_, y_), searcher)            # find aerodem points that are in neighborhood of atm points
        is_above_mindist = isempty(inds)                  # if empty -> further than mindist away from any aerodem point
        # filter values with high absolute dhdt
        return is_in_icesheet & is_above_mindist
    end
    println("Selecting flightline values...")
    filter!([:x,:y] => choose_atm, df_atm)
    # already remove some outliers here, improves standardization
    atm_to_delete = findall(abs.(df_atm.dh) .> 5 .* mad(df_atm.dh))
    deleteat!(df_atm, atm_to_delete)
    return df_atm
end

function get_aerodem_df(h_aero, bin_field_1, bin_field_2, h_ref, x, y, idx_aero)
      df_aero  = DataFrame(:x        => x[get_ix.(idx_aero, length(x))],
                         :y        => y[get_iy.(idx_aero, length(x))],
                         :h_ref    => h_ref[idx_aero],
                         :bfield_1 => bin_field_1[idx_aero],
                         :bfield_2 => bin_field_2[idx_aero],
                         :dh       => h_ref[idx_aero] - h_aero[idx_aero],
                         :h        => h_aero[idx_aero],
                         :idx      => idx_aero,
                         :source  .=> :aerodem )
    return df_aero
end

function get_itp_interp(bin_centers_1, bin_centers_2, field)
    itp  = interpolate((bin_centers_1, bin_centers_2), field, Gridded(Linear()))
    itp  = extrapolate(itp, Interpolations.Flat())
    return itp
end

function interp_nans!(bin_centers_1, field)
    b1     = get_ix.(1:length(field), length(bin_centers_1))
    b2     = get_iy.(1:length(field), length(bin_centers_1))
    gtb    = make_geotable(vec(field), b1, b2)
    out    = gtb |> InterpolateNaN(IDW())
    field .= reshape(out.Z, size(field))
    return
end

# copied and slightly edited from Geomorphometry.jl (https://github.com/Deltares/Geomorphometry.jl/blob/main/src/terrain.jl)
# not importing that package to have control over dependencies
const nbkernel = LocalFilters.Kernel{Int8,2}(reshape(1:9, 3, 3))
@inline @inbounds function horn(v, a, b)
    if b == 1
        return (v[1], v[2] + a, v[3], v[4] + a, v[5])
    elseif b == 2
        return (v[1], v[2] + 2a, v[3], v[4], v[5])
    elseif b == 3
        return (v[1], v[2] + a, v[3] + a, v[4], v[5])
    elseif b == 4
        return (v[1], v[2], v[3], v[4] + 2a, v[5])
    elseif b == 5
        return v
    elseif b == 6
        return (v[1], v[2], v[3] + 2a, v[4], v[5])
    elseif b == 7
        return (v[1] + a, v[2], v[3], v[4] + a, v[5])
    elseif b == 8
        return (v[1] + 2a, v[2], v[3], v[4], v[5])
    elseif b == 9
        return (v[1] + a, v[2], v[3] + a, v[4], v[5])
    end
end
"""
    slope(dem::Matrix{<:Real}; cellsize=1.0, method=Horn())

Slope is the rate of change between a cell and its neighbors as defined in Burrough, P. A., and McDonell, R. A., (1998, Principles of Geographical Information Systems).
"""
function slope(dem::AbstractMatrix{<:Real}; cellsize=1.0)
    dst = copy(dem)
    slope!(dst, dem, cellsize)
end
function slope!(dst, dem::AbstractMatrix{<:Real}, cellsize)  # hard-coded to Horn method
    initial(A) = (zero(eltype(A)), zero(eltype(A)), zero(eltype(A)), zero(eltype(A)), cellsize)
    store!(d, i, v) = @inbounds d[i] = atand(
        √(
            ((v[1] - v[2]) / (8 * v[5]))^2 + ((v[3] - v[4]) / (8 * v[5]))^2
        ))
    return localfilter!(dst, dem, nbkernel, initial, horn, store!)
end

bin1_fct(x, grd) = slope(x, cellsize=grd)

function get_stddization_fcts(dict_file)
    dict = load(dict_file)
    @unpack bin_centers_1, bin_centers_2, nmads, meds, std_y, mean_y = dict

    itp_var     = get_itp_interp(bin_centers_1, bin_centers_2, nmads)
    itp_bias    = get_itp_interp(bin_centers_1, bin_centers_2, meds)

    function standardize(dh, bin_field_1::Vector, bin_field_2::Vector)
        dh_detrend  = (dh .- itp_bias.(bin_field_1, bin_field_2)) ./  itp_var.(bin_field_1, bin_field_2)
        dh_detrend .= (dh_detrend .- mean_y) ./ std_y
        return dh_detrend
    end
    function destandardize(dh, bin_field_1::Vector, bin_field_2::Vector; add_mean=true)
        dh_std      = dh .* std_y .* itp_var.(bin_field_1,bin_field_2)
        if !add_mean return dh_std end
        dh_mean     = itp_bias.(bin_field_1,bin_field_2) .+ itp_var.(bin_field_1,bin_field_2) .* mean_y
        return dh_std + dh_mean
    end
    return (;standardize, destandardize)
end

function standardizing_2D(df::DataFrame; nbins1, nbins2, min_n_sample=100, fig_path)
    bin_field_1 = df[!,:bfield_1]
    bin_field_2 = df[!,:bfield_2]
    y_binned, bin_centers_1, bin_centers_2 = bin_equal_sample_size(bin_field_1, bin_field_2, Float64.(df.dh), nbins1, nbins2)
    for yb in y_binned[length.(y_binned) .< min_n_sample]
        deleteat!(yb, 1:length(yb))
        push!(yb, NaN)
    end

    # variance
    nmads       = std.(y_binned)
    Plots.heatmap(bin_centers_1, bin_centers_2, nmads')
    Plots.savefig(joinpath(fig_path, "nmads_2Dbinning_withnans.png"))
    interp_nans!(bin_centers_1, nmads)
    @assert !any(isnan.(nmads))
    itp_var     = get_itp_interp(bin_centers_1, bin_centers_2, nmads)
    # bias
    meds        = median.(y_binned)
    Plots.heatmap(bin_centers_1, bin_centers_2, meds')
    Plots.savefig(joinpath(fig_path, "medians_2Dbinning_withnans.png"))
    interp_nans!(bin_centers_1, meds)
    @assert !any(isnan.(meds))
    itp_bias    = get_itp_interp(bin_centers_1, bin_centers_2, meds)
    # standardize
    df.dh_detrend   = (df.dh .- itp_bias.(bin_field_1,bin_field_2)) ./ itp_var.(bin_field_1,bin_field_2)
    # remove outliers after standardizing
    all_to_delete = findall(abs.(df.dh_detrend) .> 10 .* mad(df.dh_detrend))
    deleteat!(df, all_to_delete)
    # make sure it's truly centered around zero and has std=1 exactly
    std_y          = std(df.dh_detrend)
    mean_y         = mean(df.dh_detrend)
    df.dh_detrend  = (df.dh_detrend .- mean_y) ./ std_y
    @printf("Kurtosis after standardization: %1.2f\n", kurtosis(df.dh_detrend))

    # return all relevant parameters so standardization function can be retrieved later
    interp_data = (;bin_centers_1, bin_centers_2, nmads, meds, std_y, mean_y)

    return df, interp_data
end

function emp_variogram(x::AbstractVector{T}, y::AbstractVector{T}, input::AbstractVector{T}; nlags=90, maxlag=7e5, fig_path="", sample_frac=1) where T <: Real
    println("Estimating variogram...")
    data = make_geotable(input, x, y)
    if sample_frac < 1
        nsamples = ceil(Int,sample_frac*length(x))
        smplr = UniformSampling(nsamples,replace=false) # replace=false ensures no point is sampled more than once
        data  = sample(data,smplr)
    end
    # compute empirical variogram
    U = data |> UniqueCoords()
    gamma = EmpiricalVariogram(U, :Z; estimator=:cressie, nlags,  maxlag) # for some reason the lsqfit has more problems fitting a good line for larger nlags
    return gamma
end

# next two functions mostly copied from https://github.com/JuliaEarth/GeoStatsFunctions.jl, modified to fit a sum of variograms with different ranges
# modified from GeoStatsFunctions.jl fit(
function fit_varg(F, f;  nVmax=3, kwargs...)
    # fit each variogram type
    res = [fit_impl(F, f, nVs; kwargs...) for nVs in 1:nVmax]
    fs, ϵs = first.(res), last.(res)

    # return best candidate
    fs[argmin(ϵs)]
end

# modified from GeoStatsFunctions.jl _fit()
function fit_impl(
    V::Type{<:Variogram},
    g::EmpiricalVariogram,
    nVs=2;
    range=nothing,
    sill=nothing,
    nugget=nothing,
    maxrange=nothing,
    maxsill=nothing,
    maxnugget=nothing
  )
    # custom ball of given radius
    ball(r) = MetricBall(r, g.distance)


    # coordinates of empirical variogram
    x = g.abscissa
    y = g.ordinate
    n = g.counts


    # discard invalid bins
    x = x[n .> 0]
    y = y[n .> 0]
    n = n[n .> 0]


    # strip units of coordinates
    ux = unit(eltype(x))
    uy = unit(eltype(y))
    x′ = ustrip.(x)
    y′ = ustrip.(y)


    # strip units of kwargs
    range′ = isnothing(range) ? range : ustrip(ux, range)
    sill′ = isnothing(sill) ? sill : ustrip(uy, sill)
    nugget′ = isnothing(nugget) ? nugget : ustrip(uy, nugget)
    maxrange′ = isnothing(maxrange) ? maxrange : ustrip(ux, maxrange)
    maxsill′ = isnothing(maxsill) ? maxsill : ustrip(uy, maxsill)
    maxnugget′ = isnothing(maxnugget) ? maxnugget : ustrip(uy, maxnugget)


    # evaluate weights
    f = nothing
    w = isnothing(f) ? n / sum(n) : map(xᵢ -> ustrip(f(xᵢ)), x)


    # objective function
    function J(θ)
      γ = sum([V(ball(rd); sill, nugget) for (rd, sill, nugget) in zip(θ[1:3:end], θ[2:3:end], θ[3:3:end])])
      sum(i -> w[i] * (γ(x′[i]) - y′[i])^2, eachindex(w, x′, y′))
    end


    # linear constraint (sill ≥ nugget)
    L(θ) = all(θ[2:3:end] .≥ θ[3:3:end]) ? 0.0 : maximum(θ[3:3:end] .- θ[2:3:end])


    # penalty for linear constraint (J + λL)
    λ = sum(yᵢ -> yᵢ^2, y′)


    # maximum range, sill and nugget
    xmax = maximum(x′)
    ymax = maximum(y′)
    rmax = isnothing(maxrange′) ? xmax : maxrange′
    smax = isnothing(maxsill′) ? ymax : maxsill′
    nmax = isnothing(maxnugget′) ? ymax : maxnugget′


    # initial guess
    rₒ = isnothing(range′) ? rmax / 3 : range′
    sₒ = isnothing(sill′) ? 0.95 * smax : sill′
    nₒ = isnothing(nugget′) ? 0.01 * smax : nugget′
    # θₒ = repeat([rₒ, sₒ, nₒ], nVs)
    θₒ = vcat([i*[rₒ, sₒ, nₒ] for i in 0.3.*(1:nVs)]...)

    # box constraints
    δ = 1e-8
    rₗ, rᵤ = isnothing(range′) ? (zero(rmax), rmax) : (range′ - δ, range′ + δ)
    sₗ, sᵤ = isnothing(sill′) ? (zero(smax), smax) : (sill′ - δ, sill′ + δ)
    nₗ, nᵤ = isnothing(nugget′) ? (zero(nmax), nmax) : (nugget′ - δ, nugget′ + δ)
    l = repeat([rₗ, sₗ, nₗ], nVs)
    u = repeat([rᵤ, sᵤ, nᵤ], nVs)


    # solve optimization problem
    sol = Optim.optimize(θ -> J(θ) + λ * L(θ), l, u, θₒ)
    ϵ = Optim.minimum(sol)
    θ = Optim.minimizer(sol)


    # optimal variogram (with units)
    # γ = V(ball(θ[1] * ux), sill=θ[2] * uy, nugget=θ[3] * uy)
    γ = sum([V(ball(rd); sill, nugget) for (rd, sill, nugget) in zip(θ[1:3:end].*ux, θ[2:3:end].*uy, θ[3:3:end].*uy)])


    γ, ϵ
end

# define variogram function to fit
function get_var(gamma; adjust_sill=true)
    varg = fit_varg(MaternVariogram, gamma, maxnugget=0.005, nVmax=2)
    vfct = typeof(varg) <: NestedVariogram ? typeof(varg.γs[1]).name.wrapper : typeof(varg).name.wrapper
    if adjust_sill
        if typeof(varg) <: NestedVariogram
            sills = [γ.sill for γ in varg.γs]
            varg = sum([vfct(γ.ball; sill=γ.sill / sum(sills) , nugget=γ.nugget) for γ in varg.γs])
        else
            varg = vfct(varg.ball; sill=1.0, nugget=varg.nugget)
        end
    end
    return varg
end

function uncertainty_from_cv(dh_binned, bin_centers, dem_ref)
    itp         = interpolate(std.(dh_binned), BSpline(Quadratic(Interpolations.Flat(OnCell()))))
    bin_c_range = range(extrema(bin_centers)..., length(bin_centers)) # convert array to range
    itp = Interpolations.scale(itp, bin_c_range)
    itp = extrapolate(itp, Interpolations.Flat())
    # save in matrix
    rec_errors          = zeros(size(dem_ref)) .+ no_data_value
    id_surf             = findall(nomissing(dem_ref .!= no_data_value, false) .|| .!ismissing.(dem_ref))
    rec_errors[id_surf] = itp.(dem_ref[id_surf])
    return itp, rec_errors
end

# for validation
function step_through_folds(flds, evaluate_fun, geotable; save_distances=false, save_coords=false)
    dif_blocks     = [Float64[] for i in flds]
    if save_distances
        NNdist_blocks  = [Float64[] for i in flds]
    end
    if save_coords
        xcoord_blocks  = [Float64[] for i in flds]
        ycoord_blocks  = [Float64[] for i in flds]
    end
    @showprogress for (j,fs) in enumerate(flds)
        # find the neighbors that the folds routine (https://github.com/JuliaEarth/GeoStatsBase.jl/blob/master/src/folding/block.jl) leaves out
        # there might be a mistake in the partitioning routine in Meshes.jl, the neighbors don't make sense (also not tested well)
        neighbors = Vector(1:length(geotable.geometry))
        deleteat!(neighbors, unique(sort([fs[1];fs[2]])))
        append!(fs[1],neighbors)      # add the removed neighbors back to training data (for Blockfolding)
        # append!(fs[2],neighbors)     # add the removed neighbors back to test data (for Ballfolding)

        sdat  = view(geotable, fs[1])
        stest = view(domain(geotable), fs[2])
        @assert length(sdat.Z) > length(stest)

        y_pred = evaluate_fun(fs[1],fs[2])
        dif_blocks[j] = mean.(y_pred) .- geotable.Z[fs[2]]

        if save_distances
            # save distance of each test point to closest training/data point
            ids_nn  = zeros(Int,length(stest))
            dist_nn = zeros(length(stest))
            for (i,st) in enumerate(stest)
                id_st, dist = searchdists(st, KNearestSearch(domain(sdat),1))
                ids_nn[i]  = id_st[1]
                dist_nn[i] = dist[1]
            end
            NNdist_blocks[j] = dist_nn
        end
        if save_coords
            # save coordinates of the test points
            crds = GeoStats.coords.(stest)
            xcoord_blocks[j] = [ustrip(c.x) for c in crds]
            ycoord_blocks[j] = [ustrip(c.y) for c in crds]
        end
    end
    rt = (vcat(dif_blocks...),)
    if save_distances
        rt = (rt..., vcat(NNdist_blocks...))
    end
    if save_coords
        rt = (rt..., vcat(xcoord_blocks...), vcat(ycoord_blocks...))
    end
    return rt
end
