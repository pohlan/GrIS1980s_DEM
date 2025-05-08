function prepare_obs(target_grid, outline_shp_file; blockspacing=400, nbins1=40, nbins2=50, coreg_grid=150)
    # define names of output directories
    main_output_dir  = joinpath("output","data_preprocessing")
    fig_path         = joinpath(main_output_dir, "figures")
    mkpath(fig_path)

    # check if files exist already
    csv_dest = joinpath(main_output_dir, "input_data_gr$(target_grid).csv")
    dict_dest = joinpath(main_output_dir, "params_gr$(target_grid).jld2")
    if all(isfile.([csv_dest, dict_dest]))
        return csv_dest, dict_dest
    end

    # get filenames at specified grid resolution
    bedmachine_original, bedm_file         = create_bedmachine_grid(target_grid)
    ref_coreg_file_ellips, ref_coreg_file_geoid, ref_file = create_grimpv2(target_grid, coreg_grid, bedmachine_original)
    aerodem_g150, obs_aero_file     = create_aerodem(target_grid, coreg_grid, outline_shp_file, bedmachine_original, ref_coreg_file_geoid)
    atm_dh_file                     = get_atm_dh_file(ref_coreg_file_ellips, ref_coreg_file_geoid, outline_shp_file, blockspacing)
    mask_file                       = create_outline_mask(target_grid, outline_shp_file, aerodem_g150)

    # read in
    h_ref        = NCDataset(ref_file)["Band1"][:,:]
    x            = NCDataset(ref_file)["x"][:]
    y            = NCDataset(ref_file)["y"][:]
    h_aero       = NCDataset(obs_aero_file)["Band1"][:,:]
    glacier_mask = NCDataset(mask_file)["Band1"][:,:]
    bedm_mask    = NCDataset(bedm_file)["mask"][:,:]
    # In areas where h_ref is NaN but h_aero is not, h_ref should be 0 instead, otherwise dh is a NaN
    id_zero = findall(ismissing.(h_ref) .&& .!ismissing.(h_aero))
    h_ref[id_zero] .= 0.0
    # don't consider the couple of negative values of h_aero
    h_aero[nomissing(h_aero,NaN) .< 0.0] .= missing

    # indices
    I_no_ocean = findall(vec(nomissing(.!ismissing.(h_ref) .&& .!ismissing.(glacier_mask) .&& (bedm_mask .!= 1.0) .&& .!(nomissing(abs.(h_ref .- h_aero),NaN)  .< 0.0), false)))
    idx_aero   = findall(vec(nomissing(.!ismissing.(h_ref) .&& .!ismissing.(glacier_mask) .&& (bedm_mask .!= 1.0) .&& .!(nomissing(abs.(h_ref .- h_aero),NaN)  .< 0.0), false) .&& .!ismissing.(h_aero)))

    # calculate bin_field on h_ref
    m_href        = nomissing(h_ref, NaN)    # bin1_fct doesn't accept type missing
    bin_field_1   = GrIS1980s_DEM.bin1_fct(m_href, target_grid)
    # calculate bin_field on h_aero to fill gaps at the terminus where h_ref is NaN (because glaciers have retreated)
    m_haero       = nomissing(h_aero, NaN)   # bin1_fct doesn't accept type missing
    b1_haero      = GrIS1980s_DEM.bin1_fct(m_haero, target_grid)
    # fill those gaps with b1_haero
    idx_aero_b1   = findall(vec(m_href .== 0.0 .&& .!isnan.(b1_haero)))
    bin_field_1[idx_aero_b1] .= b1_haero[idx_aero_b1]
    # there are still some NaN values at the edges where grimpv2 is also NaN, but h_aero has a value -> interpolate
    x_int  = x[GrIS1980s_DEM.get_ix.(I_no_ocean, length(x))]
    y_int  = y[GrIS1980s_DEM.get_iy.(I_no_ocean, length(x))]
    gtb    = GrIS1980s_DEM.make_geotable(bin_field_1[I_no_ocean], x_int, y_int)
    out    = gtb |> InterpolateNaN(IDW())
    i_nans = findall(isnan.(bin_field_1[I_no_ocean]))
    bin_field_1[I_no_ocean[i_nans]] .= out.Z[i_nans]
    # bin_field_2
    bin_field_2 = copy(m_href)
    # save as netcdf
    @assert !any(isnan.(bin_field_1[I_no_ocean]))
    @assert !any(isnan.(bin_field_2[I_no_ocean]))
    bin_field_1[isnan.(bin_field_1)] .= no_data_value
    bin_field_2[isnan.(bin_field_2)] .= no_data_value
    bfields_file  = joinpath(main_output_dir, "bin_fields_g$(target_grid).nc")
    save_netcdf(bfields_file, ref_file, [bin_field_1, bin_field_2], ["bf_1", "bf_2"], Dict("bf_1" => Dict{String,Any}(), "bf_2" => Dict{String,Any}()))
    href_file     = joinpath(main_output_dir, "h_ref_g$(target_grid).nc")
    # save h_ref as well, with the modifications
    h_ref = nomissing(h_ref, no_data_value)
    save_netcdf(href_file, ref_file, [h_ref], ["surface"], Dict("surface" => Dict{String,Any}()))

    # aerodem
    df_aero = get_aerodem_df(h_aero, bin_field_1, bin_field_2, h_ref, x, y, idx_aero)

    # atm
    df_atm  = get_ATM_df(atm_dh_file, x, y, bin_field_1, df_aero; I_no_ocean, mindist=target_grid)

    # merge aerodem and atm data
    df_all = vcat(df_aero, df_atm, cols=:union)

    # standardization
    df_all, interp_data = standardizing_2D(df_all; nbins1, nbins2, fig_path);

    # variogram
    gamma = emp_variogram(df_all.x, df_all.y, df_all.dh_detrend; maxlag=5e5, nlags=500, fig_path)

    # save
    CSV.write(csv_dest, df_all)
    to_save = (; interp_data..., gamma, I_no_ocean, idx_aero, bfields_file, href_file, coreg_grid)
    jldsave(dict_dest; to_save...)

    # make sure standardize and destandardize functions are correct
    standardize, destandardize = get_stddization_fcts(dict_dest)
    @assert all(df_all.dh_detrend .≈ standardize(df_all.dh, df_all.bfield_1, df_all.bfield_2))
    @assert all(df_all.dh         .≈ destandardize(df_all.dh_detrend, df_all.bfield_1, df_all.bfield_2))
    @assert all(df_all.h          .≈ df_all.h_ref .- destandardize(df_all.dh_detrend, df_all.bfield_1, df_all.bfield_2))
    return csv_dest, dict_dest
end

function do_kriging(output_geometry::Domain, geotable_input::AbstractGeoTable, varg::Variogram; maxn::Int)
    model  = Kriging(varg)
    interp = geotable_input |> InterpolateNeighbors(output_geometry, model, maxneighbors=maxn, prob=true)
    return interp
end

function geostats_interpolation(grd,         # make kriging a bit faster by doing it at lower resolution, then upsample back to target resolution
                                outline_shp_file,
                                csv_preprocessing, jld2_preprocessing;
                                maxn::Int)                      # maximum neighbors for interpolation method

    # define names of output directories
    main_output_dir  = joinpath("output","reconstructions")
    fig_dir          = joinpath(main_output_dir, "figures")
    mkpath(fig_dir)

    # get I_no_ocean, (de-)standardization functions and variogram from pre-processing
    df_all = CSV.read(csv_preprocessing, DataFrame)
    dict   = load(jld2_preprocessing)
    @unpack I_no_ocean, idx_aero, gamma, href_file, coreg_grid = dict
    @unpack destandardize = get_stddization_fcts(jld2_preprocessing)
    varg = get_var(gamma)

    # get filenames at grd
    bedmachine_original, _     = create_bedmachine_grid(grd)
    _, ref_coreg_file_geoid, _ = create_grimpv2(grd, coreg_grid, bedmachine_original)
    _, obs_aero_file           = create_aerodem(grd, coreg_grid, outline_shp_file, bedmachine_original, ref_coreg_file_geoid)

    # derive indices for cells to interpolate
    ir_sim      = setdiff(I_no_ocean, idx_aero)  # indices that are in I_no_ocean but not in idx_aero
    x           = NCDataset(obs_aero_file)["x"][:]
    y           = NCDataset(obs_aero_file)["y"][:]
    grid_output = PointSet([Point(xi,yi) for (xi,yi) in zip(x[get_ix.(ir_sim, length(x))], y[get_iy.(ir_sim, length(x))])])
    geotable    = make_geotable(df_all.dh_detrend, df_all.x, df_all.y)

    # prepare predicted field, fill with aerodem observations where available
    h_aero               = NCDataset(obs_aero_file)["Band1"][:,:]
    h_ref                = NCDataset(href_file)["surface"][:,:]
    h_ref                = nomissing(h_ref, 0.0)
    h_predict            = zeros(size(h_aero))
    h_predict[idx_aero] .= h_aero[idx_aero]
    bin_field_1          = bin1_fct(h_ref, grd)

    # do the kriging
    println("Kriging...")
    interp = do_kriging(grid_output, geotable, varg; maxn)
    # 'fill' aerodem with de-standardized kriging output, save as netcdf
    h_predict[ir_sim]           .= h_ref[ir_sim] .- destandardize(mean.(interp.Z), bin_field_1[ir_sim], h_ref[ir_sim])
    h_predict[h_predict .<= 0.] .= no_data_value
    # save as netcdf
    dest_file_gr_kriging         = get_rec_file_kriging(grd, maxn)
    std_uncertainty              = NCDataset(get_std_uncrt_file("kriging", grd))["std_uncertainty"][:,:]
    attributes  = Dict("surface" => Dict{String, Any}("long_name" => "ice surface elevation",
                                                      "units" => "m"),
                       "std_uncertainty" => Dict{String,Any}("long_name" => "standard deviation of error estimated from cross-validation",
                                                                 "units" => "m")
                        )
    save_netcdf(dest_file_gr_kriging, obs_aero_file, [h_predict, std_uncertainty], ["surface", "std_uncertainty"], attributes)

    # save interp.Z directly
    m_interp = zeros(size(h_predict))
    m_interp[ir_sim]            .= mean.(interp.Z)
    m_interp[m_interp .== 0.0]  .= no_data_value
    dest_m_interp                = joinpath(main_output_dir, "interpolated_dh_std_kriging.nc")
    save_netcdf(dest_m_interp, obs_aero_file, [m_interp], ["surface"], Dict("surface" => Dict{String,Any}()))
    return dest_file_gr_kriging
end
