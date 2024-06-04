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
function filter_per_bin!(y_binned; cutoff=5.0)
    println("Removing outliers per bin..")
    n_before = sum(length.(y_binned))
    is_deleted = [Int[] for i in y_binned]
    @showprogress for (i, (y_bin)) in enumerate(y_binned)
        if !isempty(y_bin)
            nmad_ = StatsBase.mad(y_bin, normalize=true)
            i_to_delete = findall(abs.(y_bin) .> cutoff*nmad_)
            sort!(unique!(i_to_delete))  # indices must be unique and sorted for keepat!
            deleteat!(y_bin,  i_to_delete)
            is_deleted[i] = i_to_delete
        end
    end
    n_after  = sum(length.(y_binned))
    perc_deleted = (n_before-n_after)*100 / n_before
    @printf "%.2f %% of points removed as outliers.\n" perc_deleted
    return is_deleted
end

function get_ATM_df(fname, dhdt, x, y, df_aero; mindist=5e4, I_no_ocean)
    df_atm = CSV.read(fname, DataFrame)
    df_atm[!,:source] .= :atm
    itp = interpolate((x, y), dhdt, Gridded(Linear()))
    df_atm[!,:dhdt] = itp.(df_atm.x, df_atm.y)
    max_dhdt = std(dhdt[dhdt .!= 0])
    function choose_atm(x_, y_, dhdt_)::Bool
        ix    = findmin(abs.(x_ .- x))[2]
        iy    = findmin(abs.(y_ .- y))[2]
        iglob = get_global_i(ix, iy, length(x))
        # filter values outside the ice sheet
        is_in_icesheet = iglob ∈ I_no_ocean
        # filter values overlapping or close to aerodem
        dist_to_aero = minimum(pairwise(Distances.euclidean, [x_ y_], [df_aero.x df_aero.y], dims=1)[:])
        is_above_mindist = dist_to_aero > mindist
        # filter values with high absolute dhdt
        has_low_dhdt = 0.0 < abs(dhdt_) < max_dhdt
        return is_in_icesheet & is_above_mindist & has_low_dhdt
    end
    println("selecting flightline values...")
    id = sort(StatsBase.sample(1:size(df_atm,1), 40000, replace=false))  # without pre-selecting a subsample of points the filtering takes a long time
    keepat!(df_atm, id)
    filter!([:x,:y,:dhdt] => choose_atm, df_atm)
    # already remove some outliers here, improves standardization
    atm_to_delete = findall(abs.(df_atm.dh) .> 5 .* mad(df_atm.dh))
    deleteat!(df_atm, atm_to_delete)
    return df_atm
end


""""
    remove_small_bins!()
Remove bins where number of samples is below a threshold.
Input
    - min_n_sample: number of samples that a bin is required to have
"""
function remove_small_bins!(A_binned::Vector...; min_n_sample=80)  # 1D
    i_rm  = findall(length.(A_binned[1]) .< min_n_sample)
    for A in A_binned
        deleteat!(A, i_rm)
    end
    @printf("Removed %d bins, %d bins left.\n", length(i_rm), length(A_binned[1]))
    return
end
function remove_small_bins(bin_centers_1::Vector, bin_centers_2::Vector, A_binned::Matrix; min_n_sample=80) # 2D
    i_rm = findall(length.(A_binned) .< min_n_sample)
    ix_rm = sort(unique([first(i_rm[i].I) for i in eachindex(i_rm)]))
    iy_rm = sort(unique([last(i_rm[i].I) for i in eachindex(i_rm)]))
    if length(ix_rm) < length(iy_rm) || (length(ix_rm) == length(iy_rm) && size(A_binned,1)>size(A_binned,2))
        idx_x = Vector(1:size(A_binned, 1))
        deleteat!(idx_x, ix_rm)
        A_binned = A_binned[idx_x, :]
        # idx_binned = idx_binned[idx_x, :]
        bin_centers_1 = bin_centers_1[idx_x]
        @printf("Removed %d rows, %d rows left.\n", length(ix_rm), length(idx_x))
    else
        idx_y = Vector(1:size(A_binned, 2))
        deleteat!(idx_y, iy_rm)
        A_binned = A_binned[:, idx_y]
        # idx_binned = idx_binned[:, idx_y]
        bin_centers_2 = bin_centers_2[idx_y]
        @printf("Removed %d columns, %d columns left.\n", length(iy_rm), length(idx_y))
    end
    @assert all(length.(A_binned) .> min_n_sample)
    return bin_centers_1, bin_centers_2, A_binned
end

function standardizing_2D(df::DataFrame, bfield1::Symbol, bfield2::Symbol; nbins1, nbins2, min_n_sample=100, fig_path="")
    bin_field_1 = df[!,bfield1]
    bin_field_2 = df[!,bfield2]
    y_binned, bin_centers_1, bin_centers_2 = bin_equal_sample_size(bin_field_1, bin_field_2, Float64.(df.dh), nbins1, nbins2)
    filter_per_bin!(y_binned)
    bin_centers_1, bin_centers_2, y_binned = remove_small_bins(bin_centers_1, bin_centers_2, y_binned; min_n_sample)
    # variance
    nmads       = mad.(y_binned, normalize=true)
    itp_mad_lin = interpolate((bin_centers_1, bin_centers_2), nmads, Gridded(Linear()))
    itp_mad_lin = extrapolate(itp_mad_lin, Interpolations.Flat())
    # bias
    meds        = median.(y_binned)
    itp_med_lin = interpolate((bin_centers_1, bin_centers_2), meds, Gridded(Linear()))
    itp_med_lin = extrapolate(itp_med_lin, Interpolations.Flat())
    # standardize
    df.dh_detrend   = (df.dh .- itp_med_lin.(bin_field_1,bin_field_2)) ./ itp_mad_lin.(bin_field_1,bin_field_2)
    # remove outliers again after standardizing
    all_to_delete = findall(abs.(df.dh_detrend) .> 4 .* mad(df.dh_detrend))
    deleteat!(df, all_to_delete)
    # make sure it's truly centered around zero and has std=1 exactly
    std_y          = std(df.dh_detrend)
    mean_y         = mean(df.dh_detrend)
    df.dh_detrend  = (df.dh_detrend .- mean_y) ./ std_y
    if !isempty(fig_path)
        # nmads
        Plots.heatmap(bin_centers_1, bin_centers_2, nmads')
        Plots.savefig(joinpath(fig_path, "nmads_2Dbinning.png"))
        # bias
        Plots.heatmap(bin_centers_1, bin_centers_2, meds')
        Plots.savefig(joinpath(fig_path, "medians_2Dbinning.png"))
        # histograms
        Plots.histogram(df.dh_detrend, label="Standardized observations", xlims=(-10,10), normalize=:pdf, nbins=1000, wsize=(600,500), linecolor=nothing)
        Plots.plot!(Normal(), lw=1, label="Normal distribution", color="black")
        Plots.savefig(joinpath(fig_path,"histogram_standardization.png"))
    end
    return df, std_y, mean_y
end

function make_geotable(input::Vector, x::Vector, y::Vector)
    table    = (;Z = input)
    coords   = [(xi, yi) for (xi,yi) in zip(x,y)]
    geotable = georef(table, coords)
    return geotable
end

function fit_variogram(x::Vector{T}, y::Vector{T}, input::Vector{T}; nlags=90, maxlag=7e5, custom_var, param_cond, p0, fig_path="", sample_frac=1) where T <: Real
    println("Estimating variogram...")
    # nx, ny = length(x), length(y)
    # X = repeat(x, 1, ny)
    # X_data = reshape(X, nx*ny, 1)[idx_outer[idx_inner]]
    # Y = repeat(y', nx, 1)
    # Y_data = reshape(Y, nx*ny, 1)[idx_outer[idx_inner]]
    # table = (Z = z[idx_inner], )
    # coords = [(X_data[i], Y_data[i]) for i in 1:length(X_data)]
    data = make_geotable(input, x, y)
    if sample_frac < 1
        nsamples = ceil(Int,sample_frac*length(x))
        smplr = UniformSampling(nsamples,replace=false) # replace=false ensures no point is sampled more than once
        data  = sample(data,smplr)
    end
    # compute empirical variogram
    U = data |> UniqueCoords()
    gamma = EmpiricalVariogram(U, :Z; estimator=:cressie, nlags,  maxlag) # for some reason the lsqfit has more problems fitting a good line for larger nlags
    # gi    = findall(gamma.ordinate .!= 0)   # sometimes there are γ=0 values, not sure why
    function get_γ(x, params)
        if !param_cond(params)  # some parameters are not accepted
            return -9999.0 .* ones(length(x))
        end
        f = custom_var(params)
        return f.(x)
    end
    # fit a covariance function
    ff = LsqFit.curve_fit(get_γ, gamma.abscissa, gamma.ordinate, p0);
    println(ff.param)
    if abs(sum(ff.param[4:6]) - 1) < 0.1 @warn "Revisit standardization, variogram doesn't seem to reach variance ~1.0." end
    varg = custom_var(ff.param)

    # plot
    if !isempty(fig_path)
        Plots.scatter(gamma.abscissa .* 1e-3, gamma.ordinate, label="Empirical variogram", color=:black, markerstrokewidth=0, wsize=(1400,800), xlabel="Distance [km]", bottommargin=10Plots.mm, leftmargin=4Plots.mm)
        Plots.plot!([0;gamma.abscissa] .* 1e-3, varg.([0;gamma.abscissa]), label="Spherical variogram fit", lw=2)
        Plots.savefig(joinpath(fig_path,"variogram.png"))
    end
    @assert varg.(1:10:100) != map(custom_var(p0),[1:10:100]...) "Fitting of variogram failed, choose better initial parameters or reduce nlags."
    return varg, ff
end


function generate_random_fields(;std_devs, corr_ls, x, y, std_z, rec, template_file, n_fields)
    k_m = 100.0
    nh  = 10000
    lx = x[end] - x[1]
    ly = y[end] - y[1]
    nx, ny = length(x), length(y)
    dest_files = Vector{String}(undef, n_fields)
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
        dest_files[i]  = "output/rec_rand_id_$(i).nc"
        svd_IceSheetDEM.save_netcdf(dest_files[i], template_file, [rftot_smooth], ["surface"], Dict("surface" => Dict()))
    end
    return dest_files
end



function residual_analysis(rec_file, bedm_file, obs_aero_file, obs_ATM_file,
                           dh_obs_long_file, dh_obs_short_file, n_years_short;
                           dh_thresh=1.0, n_bins_1=11, n_bins_2=50, n_fields=10, do_figures=false)  # 2D
    # read in
    rec          = NCDataset(rec_file)["surface"][:]
    x            = NCDataset(bedm_file)["x"][:]
    y            = NCDataset(bedm_file)["y"][:]
    mask         = NCDataset(bedm_file)["mask"][:]
    obs_aero     = NCDataset(obs_aero_file, "Band1")
    obs_ATM      = NCDataset(obs_ATM_file, "surface")
    obs_GrIMP    = NCDataset(bedm_file, "surface")
    dh_obs_short = NCDataset(dh_obs_short_file, "Band1")
    dh_obs_long  = NCDataset(dh_obs_long_file,  "Band1")

    # get indices of where to use each observational dataset
    idx_aero  = findall((rec .> 0) .&& (mask .!= 1.0) .&& (obs_aero  .> 0) .&& (abs.(rec .- obs_aero ) .> 0.0))
    # idx_GrIMP = findall((rec .> 0) .&& (mask .!= 1.0) .&& (obs_GrIMP .> 0) .&& (abs.(rec .- obs_GrIMP) .> 0.0) .&& .!bool_aero .&& (abs.(dh_obs_long) .< dh_thresh))
    idx_ATM   = findall((rec .> 0) .&& (mask .!= 1.0) .&& (obs_ATM   .> 0) .&& (abs.(rec .- obs_ATM  ) .> 0.0) .&& .!bool_aero .&& .!bool_GrIMP .&& (abs.(dh_obs_short)./n_years_short.*10 .< dh_thresh))

    # put the observations together in one matrix
    obs = copy(obs_aero)
    # obs[idx_GrIMP] = obs_GrIMP[idx_GrIMP]
    obs[idx_ATM]   = obs_ATM[idx_ATM]

    # calculate residual
    dh  = zeros(size(obs))
    dh[idx_aero]  = obs_aero[idx_aero]   .- rec[idx_aero]
    # dh[idx_GrIMP] = obs_GrIMP[idx_GrIMP] .- rec[idx_GrIMP]
    dh[idx_ATM]   = obs_ATM[idx_ATM]     .- rec[idx_ATM]

    idx_all = [idx_aero..., idx_ATM...]
    df_aero = DataFrame(:x      => x[get_ix.(idx_aero, length(x))],
                        :y      => y[get_iy.(idx_aero, length(x))],
                        :h      => obs[idx_aero],
                        :dhdt   => dhdt[idx_aero],
                        :dh     => dh[idx_aero],
                        :idx    => idx_aero,
                        :source => :aerodem)

    # df_atm = get_ATM(...)

    df_all = vcat(df_aero, df_atm, cols=:intersect)


    # index of points to be considered for statistical analysis  (has 'holes' in areas where GrIMP data not usable because too much elevation change)
    idx = findall(vec(obs) .!= 0 .&& vec(dh) .!= 0 .&& vec(dh_obs_long) .!= 0)
    # index of points where random field should be generated (entire ice sheet)
    ir_random_field = findall(obs .> 0 .|| rec .> 0)

    # # do binning over two variables
    # bin_field_1 = abs.(dh_obs_long)
    # bin_field_2 = rec
    # dh_binned,  bin_centers_1, bin_centers_2 = bin_equal_sample_size(bin_field_1[idx], bin_field_2[idx], dh[idx], n_bins_1, n_bins_2)
    # idx_binned, bin_centers_1, bin_centers_2 = bin_equal_sample_size(bin_field_1[idx], bin_field_2[idx], idx,     n_bins_1, n_bins_2)

    # # filter per bin
    # filter_per_bin!(idx_binned, dh_binned)
    # # remove rows or columns with bins that are too small
    # bin_centers_1, bin_centers_2, idx_binned, dh_binned = remove_small_bins(bin_centers_1, bin_centers_2, idx_binned, dh_binned)

    # # update dh and idx
    # dh_no_outliers = vcat(dh_binned...)
    # idx_no_outliers = vcat(idx_binned...)

    # # calculate nmad and do linear interpolation
    # nmads         = StatsBase.mad.(dh_binned, normalize=true)
    # itp           = interpolate((bin_centers_1, bin_centers_2), nmads, Gridded(Linear()))
    # linear_interp = extrapolate(itp, Interpolations.Flat())
    # σ_dh          = linear_interp.(bin_field_1[idx_no_outliers], bin_field_2[idx_no_outliers])

    # # calculate bias without bins
    # bias_dh  = median(dh_no_outliers)

    # # standardized, bias-corrected dh
    # z_dh     = (dh_no_outliers  .- bias_dh ) ./ σ_dh
    # std_z_dh = std(z_dh)
    # z_dh   ./= std_z_dh # normalize to std 1 (nmad doesn't do that)
    # @printf("Kurtosis when standardizing over bins: %1.2f\n", kurtosis(z_dh))

    # # standardization without binning, for comparison
    # dh_std = dh_no_outliers ./ std(dh_no_outliers)
    # @printf("Kurtosis when standardizing without bins: %1.2f\n", kurtosis(dh_std))

    df_all.dhdt = abs.(df_all.dhdt)
    df_all, std_dh_detrend, m_dh_detrend = svd_IceSheetDEM.standardizing_2D(df_all, :dhdt, :h_grimp; nbins1=n_bins_1, nbins2=n_bins_2, min_n_sample=50, fig_path)


    # variogram function for ParallelRandomFields
    function custom_var(x, params)
        γ1 = params[1]^2 .* (1 .-  exp.(-sqrt(2) * x./params[4]))
        γ2 = params[2]^2 .* (1 .-  exp.(-sqrt(2) * x./params[5]))
        γ3 = params[3]^2 .* (1 .-  exp.(-sqrt(2) * x./params[6]))
        return γ1 + γ2 + γ3
    end
    param_cond(params) = all(params[4:6].<1e6)
    p0 = [0.3, 0.3, 0.3, 1e3, 1e4, 1e5]
    varg, ff = fit_variogram(df_all.x, df_all.y, df_all.dh_detrend; maxlag=1.3e6, nlags=80, custom_var, param_cond, sample_frac=0.5, p0, fig_path)


    # ParallelRandomFields
    std_devs = ff.param[3]
    corr_ls  = ff.param[4]
    std_z = zeros(size(rec))
    std_z[ir_random_field] .= linear_interp.(bin_field_1[ir_random_field], bin_field_2[ir_random_field]) .* std_z_dh
    rf_files = generate_random_fields(;std_devs, corr_ls, x, y, std_z, rec, template_file=obs_aero_file, n_fields)

    if do_figures
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
