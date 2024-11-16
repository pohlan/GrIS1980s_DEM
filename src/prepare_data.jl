# Base.@kwdef mutable struct DataFiles

# Base.@kwdef struct CustomVariogram
#     params::Vector
#     custom_var::Function
#     varg::Variogram = custom_var(params)
# end

# function standardize(dh, bin_field_1, bin_field_2, itp_bias, itp_var)
#     dh_centered = dh .- itp_bias.(bin_field_1, bin_field_2) .- mean_y
#     dh_std      = itp_var.(bin_field_1, bin_field_2) .* std_y
#      return dh_centered ./ dh_std
# end
# function destandardize(dh, bin_field_1, bin_field_2, itp_bias, itp_var; add_mean=true)
#     dh_std      = dh .* std_y .* itp_var.(bin_field_1,bin_field_2)
#     if !add_mean return dh_std end
#     dh_mean     = itp_bias.(bin_field_1,bin_field_2) .+ mean_y
#     return dh_std + dh_mean
# end

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


# define variogram function to fit
function get_var(gamma; adjust_sill=true)
    varg = fit_varg(ExponentialVariogram, gamma, maxnugget=0.02, nVmax=2)
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

bin1_fct(x, grd) = slope(x, cellsize=grd)

function prepare_obs(target_grid, outline_shp_file; blockspacing=400, nbins1=40, nbins2=50, coreg_grid=150)
    # define names of output directories
    main_output_dir  = joinpath("output","data_preprocessing")
    fig_path         = joinpath(main_output_dir, "figures/")
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
    bin_field_1   = svd_IceSheetDEM.bin1_fct(m_href, target_grid)
    # calculate bin_field on h_aero to fill gaps at the terminus where h_ref is NaN (because glaciers have retreated)
    m_haero       = nomissing(h_aero, NaN)   # bin1_fct doesn't accept type missing
    b1_haero      = svd_IceSheetDEM.bin1_fct(m_haero, target_grid)
    # fill those gaps with b1_haero
    idx_aero_b1   = findall(vec(m_href .== 0.0 .&& .!isnan.(b1_haero)))
    bin_field_1[idx_aero_b1] .= b1_haero[idx_aero_b1]
    # there are still some NaN values at the edges where grimpv2 is also NaN, but h_aero has a value -> interpolate
    x_int  = x[svd_IceSheetDEM.get_ix.(I_no_ocean, length(x))]
    y_int  = y[svd_IceSheetDEM.get_iy.(I_no_ocean, length(x))]
    gtb    = svd_IceSheetDEM.make_geotable(bin_field_1[I_no_ocean], x_int, y_int)
    out    = gtb |> InterpolateNaN(IDW())
    i_nans = findall(isnan.(bin_field_1[I_no_ocean]))
    bin_field_1[I_no_ocean[i_nans]] .= out.Z[i_nans]
    # bin_field_2, also can't have NaN values so replace with h_aero values where applicable
    bin_field_2 = copy(m_href)
    # bin_field_2[.!ismissing.(h_aero)] .= h_aero[.!ismissing.(h_aero)]
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
    gamma = fit_variogram(df_all.x, df_all.y, df_all.dh_detrend; maxlag=5e5, nlags=500, fig_path)

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
