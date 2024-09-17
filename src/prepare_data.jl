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

    function standardize(dh, h_ref)
        bin_field_1 = bin1_fct(h_ref)
        bin_field_2 = h_ref
        dh_centered = dh .- itp_bias.(bin_field_1, bin_field_2) .- mean_y
        dh_std      = itp_var.(bin_field_1, bin_field_2) .* std_y
        return dh_centered ./ dh_std
    end
    function destandardize(dh, h_ref; add_mean=true)
        bin_field_1 = bin1_fct(h_ref)
        bin_field_2 = h_ref
        dh_std      = dh .* std_y .* itp_var.(bin_field_1,bin_field_2)
        if !add_mean return dh_std end
        dh_mean     = itp_bias.(bin_field_1,bin_field_2) .+ mean_y
        return dh_std + dh_mean
    end
    return (;standardize, destandardize)
end


# define variogram function to fit
custom_var(params) = SphericalVariogram(range=params[1], sill=params[4]) +
                     SphericalVariogram(range=params[2], sill=params[5]) +
                     SphericalVariogram(range=params[3], sill=params[6], nugget=params[7])

function prepare_obs(target_grid, outline_shp_file; blockspacing=target_grid, nbins1=7, nbins2=12)
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
    bedmachine_original, bedm_file  = create_bedmachine_grid(target_grid)
    reference_file_g150, ref_file   = create_grimpv2(target_grid, bedmachine_original)
    aerodem_g150, obs_aero_file     = create_aerodem(target_grid, outline_shp_file, bedmachine_original, reference_file_g150)
    atm_dh_file                     = get_atm_dh_file(ref_file, bedm_file, blockspacing)
    mask_file                       = create_outline_mask(target_grid, outline_shp_file, aerodem_g150)

    # read in
    h_ref        = NCDataset(ref_file)["Band1"][:]
    x            = NCDataset(ref_file)["x"][:]
    y            = NCDataset(ref_file)["y"][:]
    h_aero       = NCDataset(obs_aero_file)["Band1"][:]
    glacier_mask = NCDataset(mask_file)["Band1"][:]
    bedm_mask    = NCDataset(bedm_file)["mask"][:]

    # indices
    I_no_ocean = findall(vec(.!ismissing.(glacier_mask) .&& (bedm_mask .!= 1)))
    idx_aero   = findall(vec(.!ismissing.(h_aero) .&& (h_aero  .> 0) .&& .!ismissing.(h_ref) .&& (bedm_mask .!= 1.0) .&& .!ismissing.(glacier_mask) .&& (abs.(h_ref .- h_aero ) .> 0.0)))

    # aerodem
    df_aero = get_aerodem_df(h_aero, h_ref, x, y, idx_aero)

    # atm
    df_atm  = get_ATM_df(atm_dh_file, x, y, h_ref, df_aero, main_output_dir; I_no_ocean)

    # merge aerodem and atm data
    df_all = vcat(df_aero, df_atm, cols=:union)

    # plot before standardizing
    Plots.scatter(df_all.x, df_all.y, marker_z=df_all.dh, label="", markersize=2.0, markerstrokewidth=0, cmap=:RdBu, clims=(-50,50), aspect_ratio=1, xlims=(-7e5,8e5), xlabel="Easting [m]", ylabel="Northing [m]", colorbar_title="[m]", title="dh non-standardized", grid=false, wsize=(1700,1800))
    Plots.savefig(joinpath(fig_path,"data_non-standardized.png"))

    # standardization
    df_all, interp_data = standardizing_2D(df_all; nbins1, nbins2, fig_path);

    # plot after standardizing
    Plots.scatter(df_all.x, df_all.y, marker_z=df_all.dh_detrend, label="", markersize=2.0, markerstrokewidth=0, cmap=:RdBu, clims=(-4,4), aspect_ratio=1, xlims=(-7e5,8e5), xlabel="Easting [m]", ylabel="Northing [m]", colorbar_title="[m]", title="Standardized elevation difference (GrIMP - historic)", grid=false, wsize=(1700,1800))
    Plots.savefig(joinpath(fig_path,"data_standardized.png"))

    # variogram
    param_cond(params) = all(params .> 0) && all(params[4:6].<1.0) && all(params[1:3] .< 1.5e6) && params[7] .< 0.5
    # initial guess for parameters
    p0 = [6e4, 3e5, 9e5, 0.3, 0.5, 0.2, 0.25]
    _, ff = fit_variogram(F.(df_all.x), F.(df_all.y), F.(df_all.dh_detrend); maxlag=7e5, nlags=200, custom_var, param_cond, sample_frac=0.5, p0, fig_path)

    # force overall variance of variogram to be one ???
    # ff.param[4:6] .= ff.param[4:6] ./ sum(ff.param[4:6])
    display(custom_var(ff.param))

    # save
    CSV.write(csv_dest, df_all)
    to_save = (; interp_data..., params = ff.param, I_no_ocean, idx_aero)
    jldsave(dict_dest; to_save...)
    return csv_dest, dict_dest
end
