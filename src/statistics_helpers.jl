function make_geotable(Z_data, xs, ys)
    table    = (; Z=Float32.(Z_data))
    coords   = [(xi,yi) for (xi,yi) in zip(xs, ys) ]
    return georef(table,coords)
end

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


"""
    filter_per_bin

Input:
- idx_binned
-
"""
function filter_per_bin!(y_binned; cutoff=7.0)
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
    @printf("%.2f %% of points removed as outliers.\n", perc_deleted)
    return is_deleted
end

function replace_missing(A,c)
    A[ismissing.(A)] .= c
    return A
end

function get_ATM_df(fname, x, y, bin_field_1, df_aero, main_output_dir; mindist=1e4, I_no_ocean, force=false)
    dest_file = joinpath(main_output_dir, "df_atm_non-standardized.csv")
    if isfile(dest_file) && !force
        return CSV.read(dest_file, DataFrame)
    end
    df_atm = CSV.read(fname, DataFrame)
    df_atm[!,:source] .= :atm
    df_atm[!,:h]      .= df_atm.h_ref .- df_atm.dh
    sort!(unique!(x))   # necessary for interpolation
    sort!(unique!(y))
    itp = interpolate((x, y), bin_field_1, Gridded(Linear()))
    df_atm[!,:bfield_1] = itp.(df_atm.x, df_atm.y)
    df_atm[!,:bfield_2] = df_atm.h_ref
    # max_dhdt = std(dhdt[dhdt .!= 0])
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
        # has_low_dhdt = 0.0 < abs(dhdt_) < max_dhdt
        return is_in_icesheet & is_above_mindist # & has_low_dhdt
    end
    println("selecting flightline values...")
    filter!([:x,:y] => choose_atm, df_atm)
    # already remove some outliers here, improves standardization
    atm_to_delete = findall(abs.(df_atm.dh) .> 5 .* mad(df_atm.dh))
    deleteat!(df_atm, atm_to_delete)

    # save so it's faster next time
    CSV.write(dest_file, df_atm)
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

""""
    remove_small_bins!()
Remove bins where number of samples is below a threshold.
Input
    - min_n_sample: number of samples that a bin is required to have
"""
function remove_small_bins!(A_binned::Vector; min_n_sample=80)  # 1D
    i_rm  = findall(length.(A_binned[1]) .< min_n_sample)
    for A in A_binned
        deleteat!(A, i_rm)
    end
    @printf("Removed %d bins, %d bins left.\n", length(i_rm), length(A_binned[1]))
    return
end
function remove_small_bins(bin_centers_1::Vector, bin_centers_2::Vector, A_binned::Matrix; min_n_sample=80) # 2D
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
    @assert all(length.(A_binned) .>= min_n_sample)
    return bin_centers_1, bin_centers_2, A_binned
end

function get_itp_interp(bin_centers_1, bin_centers_2, field)
    itp  = interpolate((bin_centers_1, bin_centers_2), field, Gridded(Linear()))
    itp  = extrapolate(itp, Interpolations.Flat())
    return itp
end

function standardizing_2D(df::DataFrame; nbins1, nbins2, min_n_sample=30, fig_path)
    bin_field_1 = df[!,:bfield_1]
    bin_field_2 = df[!,:bfield_2]
    y_binned, bin_centers_1, bin_centers_2 = bin_equal_sample_size(bin_field_1, bin_field_2, Float64.(df.dh), nbins1, nbins2)
    bin_centers_1, bin_centers_2, y_binned = remove_small_bins(bin_centers_1, bin_centers_2, y_binned; min_n_sample)
    # variance
    nmads       = mad.(y_binned, normalize=true)
    itp_var     = get_itp_interp(bin_centers_1, bin_centers_2, nmads)
    # bias
    meds        = median.(y_binned)
    itp_bias    = get_itp_interp(bin_centers_1, bin_centers_2, meds)
    # standardize
    df.dh_detrend   = (df.dh .- itp_bias.(bin_field_1,bin_field_2)) ./ itp_var.(bin_field_1,bin_field_2)
    # remove outliers after standardizing
    all_to_delete = findall(abs.(df.dh_detrend) .> 9 .* mad(df.dh_detrend))
    deleteat!(df, all_to_delete)
    # make sure it's truly centered around zero and has std=1 exactly
    std_y          = std(df.dh_detrend)
    mean_y         = mean(df.dh_detrend)
    df.dh_detrend  = (df.dh_detrend .- mean_y) ./ std_y
    @printf("Kurtosis after standardization: %1.2f\n", kurtosis(df.dh_detrend))

    # plot and save figures
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
    Plots.heatmap(x1, x2, itp_var.(x1, x2')', xlabel="absolute elevation change over specified time period (m)", ylabel="Beucher gradient (-)", title="NMAD (-)")
    Plots.savefig(joinpath(fig_path,"nmad_interpolation.png"))
    # medians interpolation
    x1 = range(bin_centers_1[1], bin_centers_1[end], length=10000)
    x2 = range(bin_centers_2[1], bin_centers_2[end], length=1000)
    Plots.heatmap(x1, x2, itp_bias.(x1, x2')', cmap=:bwr, clims=(-20,20), xlabel="Beucher gradient (-)", ylabel="surface elevation (m)", title="median (m)")
    Plots.savefig(joinpath(fig_path,"median_interpolation.png"))

    # return all relevant parameters so standardization function can be retrieved later
    interp_data = (;bin_centers_1, bin_centers_2, nmads, meds, std_y, mean_y)

    return df, interp_data
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
    # i_nonzero = findall(gamma.ordinate .== 0)
    ff = LsqFit.curve_fit(get_γ, gamma.abscissa, gamma.ordinate, p0);
    varg = custom_var(ff.param)

    # plot
    if !isempty(fig_path)
        absc_nounits = sort(map(x -> x.val, gamma.abscissa))
        Plots.scatter(absc_nounits .* 1e-3, gamma.ordinate, label="Empirical variogram", color=:black, markerstrokewidth=0, wsize=(1400,800), xlabel="Distance [km]", bottommargin=10Plots.mm, leftmargin=4Plots.mm)
        Plots.plot!([1e-5;absc_nounits] .* 1e-3, varg.([1e-5;absc_nounits]), label="Variogram fit", lw=2, ylims=(0,1.1))
        Plots.savefig(joinpath(fig_path,"variogram.png"))
    end
    if p0 == ff.param @warn "Fitting of variogram failed, choose better initial parameters or reduce nlags." end
    return varg, ff
end

function generate_random_fields(output_dir; std_devs, corr_ls, x, y, destand, ir_random_field, rec, template_file, n_fields)
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
        rftot_destand[ir_random_field] = destand(rftot[ir_random_field], ir_random_field) ######### ToDo: second argument of destand is h_ref!!
        # smooth over a 5x5 pixel window
        rftot_smooth = mapwindow(median, rftot_destand, (5,5))
        # add the random field to the reconstruction
        rftot_smooth[ir_random_field]   .+= rec[ir_random_field]
        rftot_smooth[rftot_smooth .<= 0] .= no_data_value
        dest_files[i]  = joinpath(output_dir, "rec_rand_id_$(i).nc")
        save_netcdf(dest_files[i], template_file, [rftot_smooth], ["surface"], Dict("surface" => Dict{String,Any}()))
    end
    return dest_files
end

function SVD_random_fields(rec_file::String; nbins1::Int=10, nbins2::Int=30,  # amount of bins for 2D standardization
                           n_fields::Int=10)                                  # number of simulations
    # get filenames
    bedmachine_original, bedm_file  = create_bedmachine_grid(gr)
    aerodem_g150, obs_aero_file     = create_aerodem(gr, outline_shp_file, bedmachine_original, reference_file_g150)
    mask_file                       = create_imbie_mask(;gr, outline_shp_file, sample_path=aerodem_g150)
    dhdt_file, _                    = create_dhdt_grid(;gr, startyr=1994, endyr=2010)
    obs_ATM_file                    = get_atm_file()

    # prepare output directories
    main_output_dir  = joinpath("output","SVD_RF")
    fig_path         = joinpath(main_output_dir, "figures")
    sims_path        = joinpath(main_output_dir, "simulations")
    atm_dh_dest_file = joinpath(dirname(obs_ATM_file), "SVD_rec_minus_atm.csv")
    mkpath(fig_path)
    mkpath(sims_path)

    # define variogram function for ParallelRandomFields
    custom_var(params) =   x ->
        params[1] .* (1 .-  exp.(-sqrt(2) * x./params[3])) .+
        params[2] .* (1 .-  exp.(-sqrt(2) * x./params[4]))
    param_cond(params) = all(params .> 0.0) # conditions on parameters
    # initial guess for parameters
    p0 = [0.5, 0.5, 1e4, 4e5]

    # standardize and get variogram
    _, varg, ff, _, destand, _, I_no_ocean, _ = stddize_and_variogram(rec_file, bedm_file, obs_aero_file, obs_ATM_file, mask_file;
                                                                atm_dh_dest_file, fig_path, custom_var, param_cond, p0, nbins1, nbins2)
    @printf("Sum of variances in variogram: %.2f \n", sum(ff.param[1:2]))
    @printf("Correlation length scales: %d and %d \n", ff.param[3], ff.param[4])

    # ParallelRandomFields
    std_devs   = sqrt.(ff.param[1:2] ./ sum(ff.param[1:2]))
    corr_ls    = ff.param[3:4]
    obs        = ncread(obs_aero_file, "Band1")
    rec        = ncread(rec_file,"surface")
    x          = ncread(rec_file,"x")
    y          = ncread(rec_file,"y")
    ir_random_field = findall(vec(obs .> 0 .|| rec .> 0))
    rf_files = generate_random_fields(sims_path; std_devs, corr_ls, x, y, destand, ir_random_field, rec, template_file=obs_aero_file, n_fields)
    return rf_files
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

        # if j == 4
        #     x_ = [first(stest.domain.geoms[i[1]].coords) for i in fs[1]]
        #     y_ = [last(stest.domain.geoms[i[1]].coords)  for i in fs[1]]
        #     Plots.scatter(x_, y_, markersize=2, markerstrokewidth=0, color="grey")
        #     x_ = [first(stest.domain.geoms[i[1]].coords) for i in fs[2]]
        #     y_ = [last(stest.domain.geoms[i[1]].coords)  for i in fs[2]]
        #     Plots.scatter!(x_, y_, markersize=5, markerstrokewidth=0)
        #     Plots.savefig("foldings.png")
        #     break
        # end
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
            crds = coordinates.(stest)
            xcoord_blocks[j] = first.(crds)
            ycoord_blocks[j] = last.(crds)
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

function generate_sequential_gaussian_sim(output_geometry, geotable_input, varg; n_fields, maxn)
    method  = SEQMethod(maxneighbors=maxn)
    process = GaussianProcess(varg)
    tic   = Base.time()
    sims    = rand(process, output_geometry, geotable_input, n_fields, method)
    toc   = Base.time() - tic
    @printf("SGS took %d minutes. \n", toc / 60)
    return sims
end

function do_kriging(output_geometry::Domain, geotable_input::AbstractGeoTable, varg::Variogram; maxn::Int)
    model  = Kriging(varg)
    interp = geotable_input |> InterpolateNeighbors(output_geometry, model, maxneighbors=maxn, prob=true)
    return interp
end

function geostats_interpolation(grid_kriging, grid_out,         # make kriging a bit faster by doing it at lower resolution, then upsample back to target resolution
                                outline_shp_file,
                                csv_preprocessing, jld2_preprocessing;
                                maxn::Int,                      # maximum neighbors for interpolation method
                                method::Symbol=:kriging,        # either :kriging or :sgs
                                n_fields::Int=10)               # number of simulations in case of method=:sgs

    # get filenames at grid_kriging
    bedmachine_original, _        = create_bedmachine_grid(grid_kriging)
    reference_file_g150, ref_file = create_grimpv2(grid_kriging, bedmachine_original)
    _, obs_aero_file              = create_aerodem(grid_kriging, outline_shp_file, bedmachine_original, reference_file_g150)

    # get filename at grid_out
    # aerodem_g150, obs_aero_file_gr_out = create_aerodem(grid_out, outline_shp_file, bedmachine_original, reference_file_g150)

    # define names of output directories
    main_output_dir  = joinpath("output","geostats_interpolation")
    fig_dir          = joinpath(main_output_dir, "figures/")
    mkpath(fig_dir)

    # get I_no_ocean, (de-)standardization functions and variogram from pre-processing
    df_all = CSV.read(csv_preprocessing, DataFrame)
    dict   = load(jld2_preprocessing)
    @unpack I_no_ocean, idx_aero, params = dict
    @unpack destandardize = get_stddization_fcts(jld2_preprocessing)
    varg = custom_var(params)

    # derive indices for cells to interpolate
    ir_sim      = setdiff(I_no_ocean, idx_aero)  # indices that are in I_no_ocean but not in idx_aero
    x           = NCDataset(obs_aero_file)["x"][:]
    y           = NCDataset(obs_aero_file)["y"][:]
    grid_output = PointSet([Point(xi,yi) for (xi,yi) in zip(x[get_ix.(ir_sim, length(x))], y[get_iy.(ir_sim, length(x))])])
    geotable    = make_geotable(df_all.dh_detrend, df_all.x, df_all.y)

    # prepare predicted field, fill with aerodem observations where available
    h_aero               = NCDataset(obs_aero_file)["Band1"][:]
    h_ref                = NCDataset(ref_file)["Band1"][:]
    h_ref                = replace_missing(h_ref, 0.0)
    h_predict            = zeros(size(h_aero))
    h_predict[idx_aero] .= h_aero[idx_aero]
    std_predict          = zeros(size(h_aero))
    bin_field_1          = bin1_fct(h_ref)

    if method == :sgs  # sequential gaussian simulations
        output_path = joinpath(main_output_dir, "SEQ_simulations/")
        mkpath(output_path)
        println("Generating sequential gaussian simulations...")
        sims = generate_sequential_gaussian_sim(grid_output, geotable, varg; n_fields, maxn)
        dest_files = Vector{String}(undef, n_fields)
        for (i,s) in enumerate(sims)
            h_predict[ir_sim]           .= h_ref[ir_sim] .- destandardize(s.Z, h_ref[ir_sim])
            h_predict[h_predict .<= 0.] .= no_data_value
            dest_files[i]  = joinpath(output_path, "rec_sgs_id_$(i).nc")
            save_netcdf(dest_files[i], obs_aero_file, [h_predict], ["surface"], Dict("surface" => Dict{String,Any}()))
        end
        return dest_files
    elseif method == :kriging
        output_path = joinpath(main_output_dir, "kriging/")
        mkpath(output_path)
        # do the kriging
        println("Kriging...")
        interp = do_kriging(grid_output, geotable, varg; maxn)
        # 'fill' aerodem with de-standardized kriging output, save as netcdf
        h_predict[ir_sim]           .= h_ref[ir_sim] .- destandardize(mean.(interp.Z), bin_field_1[ir_sim], h_ref[ir_sim])
        h_predict[h_predict .<= 0.] .= no_data_value
        # field of estimated errors
        std_predict[ir_sim]         .= destandardize(std.(interp.Z), bin_field_1[ir_sim], h_ref[ir_sim], add_mean=false)
        std_predict[h_predict .== no_data_value] .= no_data_value
        # save as netcdf
        dest_file_gr_kriging         = joinpath(output_path, "rec_kriging_g$(grid_kriging).nc")
        save_netcdf(dest_file_gr_kriging, obs_aero_file, [h_predict, std_predict], ["surface", "std_error"], Dict("surface" => Dict{String,Any}(), "std_error" => Dict{String,Any}()))

        # gdalwarp to higher resolution
        # dest_file_gr_out             = joinpath(output_path, "rec_kriging_g$(grid_out).nc")
        # gdalwarp(dest_file_gr_kriging; gr=grid_out, srcnodata=string(no_data_value), dest=dest_file_gr_out)
        # # replace aerodem values with the ones warped from higher resolution, save again as netcdf
        # h_predict_gr_out             = ncread(dest_file_gr_out, "Band1")
        # h_aero_gr_out                = ncread(obs_aero_file_gr_out, "Band1")
        # ir_aero                      = findall(h_aero_gr_out .!= no_data_value)
        # h_predict_gr_out[ir_aero]   .= h_aero_gr_out[ir_aero]
        # save_netcdf(dest_file_gr_out, obs_aero_file_gr_out, [h_predict_gr_out], ["surface"], Dict("surface" => Dict{String,Any}()))

        # save interp.Z directly
        m_interp = zeros(size(h_predict))
        m_interp[ir_sim]            .= mean.(interp.Z)
        # df_aero = df_all[df_all.source .== :aerodem,:]   # also plot aerodem data that was used
        # m_interp[df_aero.idx]       .= df_aero.dh_detrend
        dest_m_interp                = joinpath(output_path, "interpolated_dh_std_kriging.nc")
        Plots.heatmap(m_interp, cmap=:bwr, clims=(-4,4))
        Plots.savefig(joinpath(fig_dir,"interp_Z.png"))
        m_interp[m_interp .== 0.0]  .= no_data_value
        save_netcdf(dest_m_interp, obs_aero_file, [m_interp], ["surface"], Dict("surface" => Dict{String,Any}()))
        return dest_file_gr_kriging
    end
end
