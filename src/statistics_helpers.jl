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
    if nx % nsampl1 < 500
        n_bins_1 -= 1
    end
    if nx % nsampl2 < 500
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


"""
    filter_per_bin

Input:
- idx_binned
-
"""
function filter_per_bin!(idx_binned, y_binned; cutoff=5.0)
    println("Removing outliers per bin..")
    n_before = sum(length.(idx_binned))
    @showprogress for (i, (i_bin, y_bin)) in enumerate(zip(idx_binned, y_binned))
        if !isempty(y_bin)
            nmad_ = StatsBase.mad(y_bin, normalize=true)
            i_to_delete = findall(abs.(y_bin) .> cutoff*nmad_)
            sort!(unique!(i_to_delete))  # indices must be unique and sorted for keepat!
            deleteat!(i_bin,  i_to_delete)
            deleteat!(y_bin,  i_to_delete)
        end
    end
    n_after  = sum(length.(idx_binned))
    perc_deleted = (n_before-n_after)*100 / n_before
    @printf "%.2f %% of points removed as outliers.\n" perc_deleted
    return
end


""""
    remove_small_bins!()
Remove bins where number of samples is below a threshold.
Input
    - min_n_sample: number of samples that a bin is required to have
"""
function remove_small_bins(bin_centers, idx_binned, y_binned; min_n_sample=80)  # 1D
    i_rm  = findall(length.(y_binned) .< min_n_sample)
    deleteat!(y_binned,   i_rm)
    deleteat!(idx_binned,  i_rm)
    deleteat!(bin_centers, i_rm)
    @printf("Removed %d bins, %d bins left.\n", length(i_rm), length(dh_binned))
    return bin_centers, idx_binned, y_binned
end
function remove_small_bins(bin_centers_1::Vector, bin_centers_2::Vector, idx_binned::Matrix{Vector{Int64}}, y_binned::Matrix{Vector{T}}; min_n_sample=80) where T <: Real# 2D
    i_rm = findall(length.(y_binned) .< min_n_sample)
    ix_rm = sort(unique([first(i_rm[i].I) for i in eachindex(i_rm)]))
    iy_rm = sort(unique([last(i_rm[i].I) for i in eachindex(i_rm)]))
    if length(ix_rm) <= length(iy_rm)
        idx_x = Vector(1:size(y_binned, 1))
        deleteat!(idx_x, ix_rm)
        y_binned = y_binned[idx_x, :]
        idx_binned = idx_binned[idx_x, :]
        bin_centers_1 = bin_centers_1[idx_x]
        @printf("Removed %d rows, %d rows left.\n", length(ix_rm), length(idx_x))
    else
        idx_y = Vector(1:size(y_binned, 2))
        deleteat!(idx_y, iy_rm)
        y_binned = y_binned[:, idx_y]
        idx_binned = idx_binned[:, idx_y]
        bin_centers_2 = bin_centers_2[idx_y]
        @printf("Removed %d columns, %d columns left.\n", length(iy_rm), length(idx_y))
    end
    @assert all(length.(y_binned) .> min_n_sample)
    return bin_centers_1, bin_centers_2, idx_binned, y_binned
end

function sample_clustered_points(idx; part_to_sample=0.1, n_per_block=3*10^3)
    n_approx      = ceil(Int,part_to_sample*length(idx)) # approximate number of points to sample
    blocks        = 1:n_per_block:length(idx)
    block_pts     = unique(rand(blocks, ceil(Int,n_approx/n_per_block)))  # sample randomly from the blocks
    idx_out       = vcat([Vector(block_pts[i]:block_pts[i]+n_per_block) for i in eachindex(block_pts)]...)  # collect all the points from the sampled blocks
    idx_rm        = findall(idx_out .> length(idx))
    deleteat!(idx_out, idx_rm)
    return idx_out
end

function get_variogram(x, y, z, idx_outer, idx_inner; nlags=2000, maxlag=2e5)
    println("Estimating variogram...")
    # make table for GeoStats.jl
    nx, ny = length(x), length(y)
    X = repeat(x, 1, ny)
    X_data = reshape(X, nx*ny, 1)[idx_outer[idx_inner]]
    Y = repeat(y', nx, 1)
    Y_data = reshape(Y, nx*ny, 1)[idx_outer[idx_inner]]
    table = (Z = z[idx_inner], )
    coords = [(X_data[i], Y_data[i]) for i in 1:length(X_data)]
    # compute empirical variogram
    data = georef(table, coords)
    U = data |> UniqueCoords()
    gamma = EmpiricalVariogram(U, :Z; estimator=:matheron, nlags,  maxlag) # for some reason the lsqfit has more problems fitting a good line for larger nlags
    gi    = findall(gamma.ordinate .!= 0)   # sometimes there are γ=0 values, not sure why
    return gamma.abscissa[gi], gamma.ordinate[gi]
end


function generate_random_fields(;std_devs, corr_ls, x, y, std_z, rec, template_file, n_fields)
    k_m = 100.0
    nh  = 10000
    lx = x[end] - x[1]
    ly = y[end] - y[1]
    nx, ny = length(x), length(y)
    dest_files = Vector{String}(undef, n_fields)
    RF_folder  = "output/SVD_RF/"
    sim_folder = joinpath(RF_folder,"simulations/")
    mkpath(sim_folder)
    for i in 1:n_fields
        rftot = zeros(nx, ny)
        for (sf, rn) in zip(std_devs, corr_ls)
            cl = (rn, rn)
            rf = generate_grf2D(lx, ly, sf, cl, k_m, nh, nx, ny, cov_typ="expon", do_reset=false, do_viz=false);
            rftot .+= Array(rf)
        end
        # smooth over a 5x5 pixel window
        rftot_smooth = mapwindow(median, rftot, (7,7))
        # multiply with nmad and sigmas to get back variability w.r.t. binning variables
        rftot_smooth .*= std_z
        # add the random field to the reconstruction
        rftot_smooth .+= rec
        dest_files[i]  = joinpath.(sim_folder, "test_rec_rand_id_$(i).nc")
        svd_IceSheetDEM.save_netcdf(dest_files[i], template_file, [rftot_smooth], ["surface"], Dict("surface" => Dict{String,Any}()))
    end
    return RF_folder, dest_files
end



function residual_analysis(rec_file, bedm_file, obs_aero_file, obs_ATM_file,
                           dh_obs_long_file, dh_obs_short_file, n_years_short;
                           dh_thresh=1.0, n_bins_1=11, n_bins_2=50, n_fields=10, do_figures=false)  # 2D
    # read in
    rec          = ncread(rec_file, "surface")
    x            = ncread(bedm_file, "x")
    y            = ncread(bedm_file, "y")
    mask         = ncread(bedm_file, "mask")
    obs_aero     = ncread(obs_aero_file, "Band1")
    obs_ATM      = ncread(obs_ATM_file, "surface")
    obs_GrIMP    = ncread(bedm_file, "surface")
    dh_obs_short = ncread(dh_obs_short_file, "Band1")
    dh_obs_long  = ncread(dh_obs_long_file,  "Band1")

    # get indices of where to use each observational dataset
    bool_aero  = (rec .> 0) .&& (mask .!= 1.0) .&& (obs_aero  .> 0) .&& (abs.(rec .- obs_aero ) .> 0.0)
    bool_GrIMP = (rec .> 0) .&& (mask .!= 1.0) .&& (obs_GrIMP .> 0) .&& (abs.(rec .- obs_GrIMP) .> 0.0) .&& .!bool_aero .&& (abs.(dh_obs_long) .< dh_thresh)
    bool_ATM   = (rec .> 0) .&& (mask .!= 1.0) .&& (obs_ATM   .> 0) .&& (abs.(rec .- obs_ATM  ) .> 0.0) .&& .!bool_aero .&& .!bool_GrIMP .&& (abs.(dh_obs_short)./n_years_short.*10 .< dh_thresh)
    idx_aero  = findall(bool_aero)
    idx_GrIMP = findall(bool_GrIMP)
    idx_ATM   = findall(bool_ATM)

    # put the observations together in one matrix
    obs = copy(obs_aero)
    obs[idx_GrIMP] = obs_GrIMP[idx_GrIMP]
    obs[idx_ATM]   = obs_ATM[idx_ATM]

    # calculate residual
    dh  = zeros(size(obs))
    dh[idx_aero]  = obs_aero[idx_aero]   .- rec[idx_aero]
    dh[idx_GrIMP] = obs_GrIMP[idx_GrIMP] .- rec[idx_GrIMP]
    dh[idx_ATM]   = obs_ATM[idx_ATM]     .- rec[idx_ATM]

    # index of points to be considered for statistical analysis  (has 'holes' in areas where GrIMP data not usable because too much elevation change)
    idx = findall(vec(obs) .!= 0 .&& vec(dh) .!= 0 .&& vec(dh_obs_long) .!= 0)
    # index of points where random field should be generated (entire ice sheet)
    ir_random_field = findall(obs .> 0 .|| rec .> 0)

    # do binning over two variables
    bin_field_1 = abs.(dh_obs_long)
    bin_field_2 = rec
    dh_binned,  bin_centers_1, bin_centers_2 = bin_equal_sample_size(bin_field_1[idx], bin_field_2[idx], dh[idx], n_bins_1, n_bins_2)
    idx_binned, bin_centers_1, bin_centers_2 = bin_equal_sample_size(bin_field_1[idx], bin_field_2[idx], idx,     n_bins_1, n_bins_2)

    # filter per bin
    filter_per_bin!(idx_binned, dh_binned)
    # remove rows or columns with bins that are too small
    bin_centers_1, bin_centers_2, idx_binned, dh_binned = remove_small_bins(bin_centers_1, bin_centers_2, idx_binned, dh_binned)

    # update dh and idx
    dh_no_outliers = vcat(dh_binned...)
    idx_no_outliers = vcat(idx_binned...)

    # calculate nmad and do linear interpolation
    nmads         = StatsBase.mad.(dh_binned, normalize=true)
    itp           = interpolate((bin_centers_1, bin_centers_2), nmads, Gridded(Linear()))
    linear_interp = extrapolate(itp, Interpolations.Flat())
    σ_dh          = linear_interp.(bin_field_1[idx_no_outliers], bin_field_2[idx_no_outliers])

    # calculate bias without bins
    bias_dh  = median(dh_no_outliers)

    # standardized, bias-corrected dh
    z_dh     = (dh_no_outliers  .- bias_dh ) ./ σ_dh
    std_z_dh = std(z_dh)
    z_dh   ./= std_z_dh # normalize to std 1 (nmad doesn't do that)
    @printf("Kurtosis when standardizing over bins: %1.2f\n", kurtosis(z_dh))

    # standardization without binning, for comparison
    dh_std = dh_no_outliers ./ std(dh_no_outliers)
    @printf("Kurtosis when standardizing without bins: %1.2f\n", kurtosis(dh_std))

    # sample points for variogram (takes too long otherwise), clustered in a certain "block" since completely random sampling will not cover range of distances well
    idx_variogram = sample_clustered_points(idx_no_outliers, part_to_sample=0.1)

    # estimate variogram
    distances, γs = get_variogram(x, y, z_dh, idx_no_outliers, idx_variogram)

    # fit the function for ParallelRandomFields manually
    function custom_var(x, params)
        γ1 = params[1]^2 .* (1 .-  exp.(-sqrt(2) * x./params[2]))
        γ2 = params[3]^2 .* (1 .-  exp.(-sqrt(2) * x./params[4]))
        γ3 = params[5]^2 .* (1 .-  exp.(-sqrt(2) * x./params[6]))
        return γ1 + γ2 + γ3
    end
    ff = LsqFit.curve_fit(custom_var, distances, γs, [0.3, 1e3, 0.3, 1e4, 0.3,  1e5]);

    # ParallelRandomFields
    std_devs = ff.param[3]
    corr_ls  = ff.param[4]
    std_z = zeros(size(rec))
    std_z[ir_random_field] .= linear_interp.(bin_field_1[ir_random_field], bin_field_2[ir_random_field]) .* std_z_dh
    RF_folder, rf_files = generate_random_fields(;std_devs, corr_ls, x, y, std_z, rec, template_file=obs_aero_file, n_fields)

    if do_figures
        fig_dir = joinpath(RF_folder, "figures/")
        mkpath(fig_dir)
        # nmad interpolation
        x1 = range(bin_centers_1[1], bin_centers_1[end], length=10000)
        x2 = range(bin_centers_2[1], bin_centers_2[end], length=1000)
        Plots.heatmap(x1, x2, linear_interp.(x1, x2')', xaxis=:log, xlabel="absolute elevation change over specified time period (m)", ylabel="surface elevation (m)", title="NMAD (-)")
        Plots.savefig(joinpath(fig_dir,"nmad_interpolation.png"))
        # qqplot
        Plots.plot(
            StatsPlots.qqplot(StatsPlots.Normal(), z_dh, title="standardized with binning", ylims=(-8,8)),
            StatsPlots.qqplot(StatsPlots.Normal(), dh_std, title="standardized without binning", ylims=(-8,8))
        )
        Plots.savefig(joinpath(fig_dir,"qqplot.png"))
        # histogram
        p1 = Plots.histogram(z_dh, xlims=(-4,4), ylims=(0, 0.7), label="standardized residual", bins=80, normalize=true, alpha=0.5, title="With binning")
        p1 = StatsPlots.plot!(Normal(0,1), label="Normal distribution", xlims=(-4,4), ylims=(0, 0.7), legend=false)
        p2 = Plots.histogram(dh_std.-median(dh_std), xlims=(-4,4), ylims=(0, 0.7), label="standardized residual", bins=120, normalize=true, alpha=0.5, title="Without binning")
        p2 = StatsPlots.plot!(Normal(0,1), label="Normal distribution", legend=false)
        Plots.plot(p1, p2)
        Plots.savefig(joinpath(fig_dir,"histograms.png"))
        # variogram
        Plots.scatter(distances, γs, xaxis=:log, label="")
        Plots.plot!(distances, custom_var(distances, ff.param), label="LsqFit", legend=:topleft)
        Plots.savefig(joinpath(fig_dir, "variogram.png"))
    end
    return rf_files
end
