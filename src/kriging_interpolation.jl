function prepare_obs(target_grid, outline_shp_file; blockspacing=400, nbins1=5, nbins2=18, coreg_grid=150, r_aero_varg=0.2)
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
    dh_error_file                   = get_atm_aero_error_file(aerodem_g150, atm_dh_file)

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

    # variogram (use only a ratio --r_aero_varg-- of the aerodme data) otherwise variogram is dominated by aerodem
    i_aero = findall(df_all.source .== :aerodem)
    i_atm  = findall(df_all.source .== :atm)
    n_aero = round(Int, r_aero_varg * length(i_aero))
    i_sampl = [rand(1:length(i_aero), n_aero); i_atm]
    gamma = emp_variogram(df_all.x[i_sampl], df_all.y[i_sampl], df_all.dh_detrend[i_sampl]; maxlag=2e5, nlags=200, fig_path)

    # error variogram
    df_error    = get_error_df(dh_error_file, x, y, I_no_ocean)
    U           = make_geotable(df_error.dh, df_error.x, df_error.y)
    gamma_error = EmpiricalVariogram(U, :Z; estimator=:cressie, nlags=20,  maxlag=1.2e4)

    # save
    CSV.write(csv_dest, df_all)
    to_save = (; interp_data..., gamma, gamma_error, I_no_ocean, idx_aero, bfields_file, href_file, coreg_grid)
    jldsave(dict_dest; to_save...)

    # make sure standardize and destandardize functions are correct
    standardize, destandardize = get_stddization_fcts(dict_dest)
    @assert all(df_all.dh_detrend .≈ standardize(df_all.dh, df_all.bfield_1, df_all.bfield_2))
    @assert all(df_all.dh         .≈ destandardize(df_all.dh_detrend, df_all.bfield_1, df_all.bfield_2))
    @assert all(df_all.h          .≈ df_all.h_ref .- destandardize(df_all.dh_detrend, df_all.bfield_1, df_all.bfield_2))
    return csv_dest, dict_dest
end

struct CustomKernel <: AbstractGPs.Kernel
    kernel::AbstractGPs.Kernel
    r::Real
end
function (k::CustomKernel)(x, y)
    val = k.kernel(x/k.r, y/k.r)  # divide by range before applying kernel
    return val
end
function varg_to_kernel(varg::NestedVariogram)
    γ1, γ2 = varg.γs
    base_kernel = ExponentialKernel()
    kernel = F(sill(γ1)) * CustomKernel(base_kernel, range(γ1).val/3) +
             F(sill(γ2)) * CustomKernel(base_kernel, range(γ2).val/3)
    return kernel
end
function varg_to_kernel(varg::ExponentialVariogram)
    base_kernel = ExponentialKernel()
    kernel = F(sill(varg)) * CustomKernel(base_kernel, range(varg).val/3)
    return kernel
end
function do_GP(coords_input, Z_input, coords_output, kernel_signal, kernel_error, σ_obs; var=true)
    # set kernel
    f_GP = GP(ZeroMean{F}(), kernel_signal)
    # do GP
    Σ_obs        = kernelmatrix(kernel_error, coords_input) .* (σ_obs * σ_obs')
    f_GP_input   = f_GP(coords_input, Σ_obs)
    posterior_gp = posterior(f_GP_input, F.(Z_input))
    pred_fct     = var ? mean_and_var : mean
    out          = pred_fct(posterior_gp(coords_output))
    return out
end

function add_sigma_obs!(df_all, standardize=nothing)
    σ_aero, σ_atm = 10.0, 0.5
    σ_obs = zeros(F, length(df_all.x))
    σ_obs .= σ_atm; σ_obs[df_all.source .== "aerodem"] .= σ_aero
    df_all.sigma_obs = σ_obs
    if !isnothing(standardize)
        df_all.sigma_obs_detrend = standardize(σ_obs, df_all.bfield_1, df_all.bfield_2)
    end
    return
end

function geostats_interpolation(grd,
                                outline_shp_file,
                                csv_preprocessing, jld2_preprocessing;
                                ℓ_block::Real=1.5e5,  # size of block to be interpolated at once
                                δl=1.5e5)             # how much data is used outside of block (width of stripe on each side)

    # define names of output directories
    main_output_dir  = joinpath("output","reconstructions")
    fig_dir          = joinpath(main_output_dir, "figures")
    mkpath(fig_dir)

    # get I_no_ocean, (de-)standardization functions and variogram from pre-processing
    df_all     = CSV.read(csv_preprocessing, DataFrame)
    coords_obs = [F.([x,y]) for (x,y) in zip(df_all.x, df_all.y)]
    dict       = load(jld2_preprocessing)
    @unpack I_no_ocean, idx_aero, gamma, gamma_error, href_file, coreg_grid = dict
    @unpack standardize, destandardize = get_stddization_fcts(jld2_preprocessing)
    varg = get_var(gamma)
    varg_error = get_var(gamma_error, nVmax=1)
    kernel_signal = varg_to_kernel(varg)
    kernel_error  = varg_to_kernel(varg_error)

    # get filenames at grd
    bedmachine_original, _     = create_bedmachine_grid(grd)
    _, ref_coreg_file_geoid, _ = create_grimpv2(grd, coreg_grid, bedmachine_original)
    _, obs_aero_file           = create_aerodem(grd, coreg_grid, outline_shp_file, bedmachine_original, ref_coreg_file_geoid)

    # derive indices and coordinates of cells to interpolate
    ir_sim      = setdiff(I_no_ocean, idx_aero)  # indices that are in I_no_ocean but not in idx_aero
    x           = NCDataset(obs_aero_file)["x"][:]; nx = length(x)
    y           = NCDataset(obs_aero_file)["y"][:]
    x_sim       = x[get_ix.(ir_sim, nx)]
    y_sim       = y[get_iy.(ir_sim, nx)]
    coords_sim  = [F.([x,y]) for (x,y) in zip(x_sim, y_sim)]

    # divide into blocks (too much memory to do all at once)
    geotable    = georef(nothing, coords_sim)
    blocks      = partition(geotable, BlockPartition(ℓ_block))
    idcs        = indices(blocks)

    # add measurement uncertainty for AeroDEM / ATM
    add_sigma_obs!(df_all, standardize)

    # prepare predicted field, fill with aerodem observations where available
    h_aero               = NCDataset(obs_aero_file)["Band1"][:,:]
    h_ref                = NCDataset(href_file)["surface"][:,:]
    h_ref                = nomissing(h_ref, 0.0)
    h_predict            = zeros(size(h_aero))
    h_predict[idx_aero] .= h_aero[idx_aero]
    σ_predict            = zeros(size(h_aero))
    σ_predict[idx_aero] .= 10.0
    bin_field_1          = bin1_fct(h_ref, grd)
    dh_predict           = zeros(size(h_aero))
    std_predict          = zeros(size(h_aero))

    # loop through blocks to do GP
    println("Gaussian Process interpolation...")
    println("Number of blocks: $(length(idcs))")
    @showprogress Threads.@threads for i_block in idcs
        # subset of coordinates to interpolate in this block
        coords_output = coords_sim[i_block]

        # find bounding values for x and y
        x_min, x_max = extrema(first.(coords_output))
        y_min, y_max = extrema(last.(coords_output))

        # indices in terms of df_all, input data (add stripe of width δl around interpolation domain)
        i_dat = findall(x_min-δl .<= df_all.x .<= x_max+δl .&& y_min-δl .<= df_all.y .<= y_max+δl)

        # do GP and save output
        m_pred, std_pred = do_GP(coords_obs[i_dat], F.(df_all.dh_detrend[i_dat]), coords_output, kernel_signal, kernel_error, df_all.sigma_obs_detrend[i_dat])
        dh_predict[ir_sim[i_block]]  .= m_pred
        std_predict[ir_sim[i_block]] .= sqrt.(std_pred) # get std from variance
    end

    # 'fill' aerodem with de-standardized kriging output, save as netcdf
    h_predict[ir_sim]           .= h_ref[ir_sim] .- destandardize(dh_predict[ir_sim], bin_field_1[ir_sim], h_ref[ir_sim])
    h_predict[h_predict .<= 0.] .= no_data_value
    σ_predict[ir_sim]           .= destandardize(std_predict[ir_sim], bin_field_1[ir_sim], h_ref[ir_sim], add_mean=false)
    σ_predict[σ_predict .<= 0.] .= no_data_value

    # save as netcdf
    dest_file_GP = get_rec_file_GP(grd)
    attributes   = Dict("surface" => Dict{String, Any}("long_name" => "ice surface elevation",
                                                       "units" => "m"),
                        "std_uncertainty" => Dict{String,Any}("long_name" => "uncertainty estimation ...", # TODO
                                                                  "units" => "m")
                        )
    save_netcdf(dest_file_GP, obs_aero_file, [h_predict, σ_predict], ["surface", "std_uncertainty"], attributes)

    # save interpolated non-standardized fields directly
    dh_predict[dh_predict .== 0.0]  .= no_data_value
    std_predict[std_predict .== 0.] .= no_data_value
    dest_dh_predict = joinpath(main_output_dir, "interpolated_dh_std_GP.nc")
    save_netcdf(dest_dh_predict, obs_aero_file, [dh_predict, std_predict], ["dh_std", "sigma_std"], Dict("dh_std" => Dict{String,Any}(),"sigma_std" => Dict{String,Any}()))

    return dest_file_GP
end
