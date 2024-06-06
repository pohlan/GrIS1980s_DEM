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

function replace_missing(A,c)
    A[ismissing.(A)] .= c
    return A
end

function get_ATM_df(fname, dhdt0, x, y, df_aero; mindist=5e4, I_no_ocean)
    df_atm = CSV.read(fname, DataFrame)
    df_atm[!,:source] .= :atm
    dhdt = replace_missing(dhdt0, 0)
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

function get_aerodem_df(aero, ref, dhdt0, x, y, idx_aero)
    dhdt = replace_missing(dhdt0, 0)
    df_aero  = DataFrame(:x       => x[get_ix.(idx_aero, length(x))],
                         :y       => y[get_iy.(idx_aero, length(x))],
                         :h       => ref[idx_aero],
                         :dhdt    => dhdt[idx_aero],
                         :dh      => ref[idx_aero] - aero[idx_aero],
                         :idx     => idx_aero,
                         :source .=> :aerodem )
    return df_aero
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
    display(length.(A_binned))
    while any(length.(A_binned) .< min_n_sample)
        i_rm               = findall(length.(A_binned) .< min_n_sample)
        ix_rm              = [first(i_rm[i].I) for i in eachindex(i_rm)]
        iy_rm              = [last(i_rm[i].I) for i in eachindex(i_rm)]
        x_maxcount, ix_max = findmax(countmap(ix_rm))
        y_maxcount, iy_max = findmax(countmap(iy_rm))
        if x_maxcount > y_maxcount || (x_maxcount == y_maxcount && size(A_binned,1)>size(A_binned,2))
            idx_x = Vector(1:size(A_binned, 1))
            deleteat!(idx_x, ix_max)
            A_binned = A_binned[idx_x, :]
            bin_centers_1 = bin_centers_1[idx_x]
            @printf("Removed one row, %d rows left.\n", length(idx_x))
        else
            idx_y = Vector(1:size(A_binned, 2))
            deleteat!(idx_y, iy_max)
            A_binned = A_binned[:, idx_y]
            bin_centers_2 = bin_centers_2[idx_y]
            @printf("Removed one column, %d columns left.\n", length(idx_y))
        end
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
        # qqplot
        Plots.plot(
            StatsPlots.qqplot(StatsPlots.Normal(), df.dh_detrend, title="standardized with binning", ylims=(-8,8)),
            StatsPlots.qqplot(StatsPlots.Normal(), (df.dh .- mean(df.dh))./std(df.dh), title="standardized without binning", ylims=(-8,8))
        )
        Plots.savefig(joinpath(fig_path,"qqplot.png"))
        # nmad interpolation
        x1 = range(bin_centers_1[1], bin_centers_1[end], length=10000)
        x2 = range(bin_centers_2[1], bin_centers_2[end], length=1000)
        Plots.heatmap(x1, x2, itp_mad_lin.(x1, x2')', xaxis=:log, xlabel="absolute elevation change over specified time period (m)", ylabel="surface elevation (m)", title="NMAD (-)")
        Plots.savefig(joinpath(fig_path,"nmad_interpolation.png"))
    end
    @printf("Kurtosis after standardization: %1.2f\n", kurtosis(df.dh_detrend))
    destandardize(dh, bin_field_1, bin_field_2) = dh .* std_y .* itp_mad_lin.(bin_field_1,bin_field_2) .+ itp_med_lin.(bin_field_1,bin_field_2) .+ mean_y
    return df, destandardize
end

function make_geotable(input::Vector, x::Vector, y::Vector)
    table    = (;Z = input)
    coords   = [(xi, yi) for (xi,yi) in zip(x,y)]
    geotable = georef(table, coords)
    return geotable
end

function fit_variogram(x::Vector{T}, y::Vector{T}, input::Vector{T}; nlags=90, maxlag=7e5, custom_var, param_cond, p0, fig_path="", sample_frac=1) where T <: Real
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
    varg = custom_var(ff.param)
    # plot
    if !isempty(fig_path)
        Plots.scatter(gamma.abscissa .* 1e-3, gamma.ordinate, label="Empirical variogram", color=:black, markerstrokewidth=0, wsize=(1400,800), xlabel="Distance [km]", bottommargin=10Plots.mm, leftmargin=4Plots.mm)
        Plots.plot!([0;gamma.abscissa] .* 1e-3, varg.([0;gamma.abscissa]), label="Variogram fit", lw=2)
        Plots.savefig(joinpath(fig_path,"variogram.png"))
    end
    if p0 == ff.param @warn "Fitting of variogram failed, choose better initial parameters or reduce nlags." end
    return varg, ff
end

function generate_random_fields(;std_devs, corr_ls, x, y, destand, ir_random_field, rec, template_file, n_fields)
    output_dir = "output/SVD_RF/simulations/"
    mkpath(output_dir)
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
        # multiply with nmad and sigmas to get back variability w.r.t. binning variables
        rftot_destand = zeros(Float32, nx, ny)
        rftot_destand[ir_random_field] = destand(rftot, ir_random_field)
        # smooth over a 5x5 pixel window
        rftot_smooth = mapwindow(median, rftot_destand, (5,5))
        # add the random field to the reconstruction
        rftot_smooth[ir_random_field] .+= rec[ir_random_field]
        rftot_smooth[rftot_smooth .<= 0] .= no_data_value
        dest_files[i]  = joinpath(output_dir, "rec_rand_id_$(i).nc")
        svd_IceSheetDEM.save_netcdf(dest_files[i], template_file, [rftot_smooth], ["surface"], Dict("surface" => Dict()))
    end
    return dest_files
end


function prepare_random_sims(ref_file, bedm_file, obs_aero_file, obs_ATM_file, dhdt_file, mask_file;
                             dh_atm_file, fig_path, custom_var, param_cond, p0, nbins1, nbins2, min_n_sample=100)
    # read in
    ref          = NCDataset(ref_file)["surface"][:]
    x            = NCDataset(ref_file)["x"][:]
    y            = NCDataset(ref_file)["y"][:]
    h_aero       = NCDataset(obs_aero_file)["Band1"][:]
    dhdt         = NCDataset(dhdt_file)["Band1"][:]
    glacier_mask = NCDataset(mask_file)["Band1"][:]
    bedm_mask    = NCDataset(bedm_file)["mask"][:]

    # indices
    I_no_ocean = findall(vec(.!ismissing.(glacier_mask) .&& (bedm_mask .!= 1)))
    idx_aero   = findall(vec(.!ismissing.(h_aero) .&& .!ismissing.(ref) .&& (bedm_mask .!= 1.0) .&& (h_aero  .> 0) .&& (abs.(ref .- h_aero ) .> 0.0)))

    # aerodem
    df_aero = get_aerodem_df(h_aero, ref, dhdt, x, y, idx_aero)

    # atm
    ref_surface = get_surface_file(ref_file, bedm_file, remove_geoid=true)
    # interpolate GrIMP DEM on ATM points and calculate difference
    py_point_interp(ref_surface, obs_ATM_file, dh_atm_file)
    # read in and select values
    df_atm = get_ATM_df(dh_atm_file, dhdt, x, y, df_aero; I_no_ocean)

    # merge aerodem and atm data
    df_all = vcat(df_aero, df_atm, cols=:intersect)

    # standardize, describing variance and bias as a function of dhdt and elevation
    df_all.dhdt = abs.(df_all.dhdt)
    df_all, destandardize = standardizing_2D(df_all, :dhdt, :h; nbins1, nbins2, min_n_sample, fig_path);
    dhdt = replace_missing(dhdt, 0.0)   # the variance and bias linear interpolation function don't accept a type missing
    ref  = replace_missing(ref, 0.0)
    destand(dh, idx) = destandardize(dh[idx], abs.(dhdt[idx]), ref[idx])

    # plot after standardizing
    Plots.scatter(df_all.x, df_all.y, marker_z=df_all.dh_detrend, label="", markersize=2.0, markerstrokewidth=0, cmap=:RdBu, clims=(-4,4), aspect_ratio=1, xlims=(-7e5,8e5), xlabel="Easting [m]", ylabel="Northing [m]", colorbar_title="[m]", title="Standardized elevation difference (GrIMP - historic)", grid=false, wsize=(1700,1800))
    Plots.savefig(joinpath(fig_path,"data_standardized.png"))

    # variogram
    varg, ff = fit_variogram(df_all.x, df_all.y, df_all.dh_detrend; maxlag=1.2e6, nlags=200, custom_var, param_cond, sample_frac=0.05, p0, fig_path)
    return df_all, varg, ff, destand, I_no_ocean, idx_aero
end

function SVD_random_fields(rec_file, bedm_file, obs_aero_file, obs_ATM_file, dhdt_file, mask_file;
                           nbins1=10, nbins2=30, n_fields=10)
    fig_path = "output/SVD_RF/figures/"
    mkpath(fig_path)
    dh_atm_file = joinpath(dirname(obs_ATM_file), "SVD_rec_minus_atm.csv")

    # define variogram function for ParallelRandomFields
    custom_var(params) =   x ->
        params[1]^2 .* (1 .-  exp.(-sqrt(2) * x./params[3])) .+
        params[2]^2 .* (1 .-  exp.(-sqrt(2) * x./params[4]))
    param_cond(params) = all(0.0 .< params[1:2] .< 1.0) && all(sum(params[1:2].^2).>0.95) && all(params[3:4] .< 1e6) # conditions on parameters
    p0 = [0.8, 0.5, 1e4, 4e5]

    # standardize and get variogram
    _, _, ff, destand, I_no_ocean, _ = prepare_random_sims(rec_file, bedm_file, obs_aero_file, obs_ATM_file, dhdt_file, mask_file;
                                                    dh_atm_file, fig_path, custom_var, param_cond, p0, nbins1, nbins2)

    # ParallelRandomFields
    std_devs   = ff.param[1:2]
    corr_ls    = ff.param[3:4]
    obs        = ncread(obs_aero_file, "Band1")
    rec        = ncread(rec_file,"surface")
    x          = ncread(rec_file,"x")
    y          = ncread(rec_file,"y")
    ir_random_field = findall(vec(obs .> 0 .|| rec .> 0))
    rf_files = generate_random_fields(;std_devs, corr_ls, x, y, destand, ir_random_field, rec, template_file=obs_aero_file, n_fields)
    return rf_files
end
