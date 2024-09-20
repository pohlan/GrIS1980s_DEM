function read_model_data(;which_files=nothing,       # indices of files used for training     ; e.g. 1:10, default all available
                          tsteps=nothing,            # indices of time steps used for training; e.g.  ""      ""
                          model_files,
                          I_no_ocean)
    println("Reading in model data...")
    # determine indices for files
    if which_files === nothing
        which_files = 1:length(model_files)
    elseif which_files[end] > length(model_files)
        error("Number of files out of bound.")
    end
    files_out  = model_files[which_files]
    nf         = length(files_out)

    # determine total number of time steps
    nts = []
    for f in files_out
        ds   = NCDataset(f)
        nt_f = size(ds["usurf"], 3)
        push!(nts, nt_f)
        close(ds)
    end
    if isnothing(tsteps)
        nttot = sum(nts)
    else
        nttot = sum(min.(nts, length(tsteps)))
    end

    # determine number of cells in x and y direction
    ds = NCDataset(files_out[1])
    nx, ny = size(ds["usurf"])[1:2]
    close(ds)
    # build data matrix
    Data = zeros(F, length(I_no_ocean), nttot)
    ntcount = 0
    @showprogress for (k, file) in enumerate(files_out)
    d = ncread(file, "usurf")
    if isnothing(tsteps)
        ts = 1:size(d,3)
    elseif minimum(tsteps) > size(d,3)
        continue
    else
        ts = tsteps[1]:min(tsteps[end], size(d, 3))
    end
    nt_out = length(ts)
    data = reshape(d[:,:,ts], ny*nx, nt_out)
    @views Data[:, ntcount+1 : ntcount+nt_out] = data[I_no_ocean,:]
    ntcount += nt_out
    end
    return Data
end

function prepare_model(model_files, standardize, h_ref, I_no_ocean, r, output_dir; use_arpack=false, input="dh_detrend")
    # load model data and calculate difference to reference DEM
    Data_ice  = read_model_data(;model_files,I_no_ocean)

    # calculate difference to reference DEM
    data_ref  = h_ref[I_no_ocean]
    if input !== "h"
        Data_ice .= data_ref .- Data_ice
    end

    # standardize model data
    if input == "dh_detrend"
        for dc in eachcol(Data_ice)
            dc .= standardize(dc, data_ref)
        end
    end
    # subtract mean again for it to be truly centered
    data_mean  = mean(Data_ice, dims=2)
    Data_ice  .= Data_ice .- data_mean

    # compute SVD
    if use_arpack
        nsv = min(r, size(Data_ice, 2)-1) # actual truncation is later, but takes too long if r is unnecessarily high here
        B, _ = svds(Data_ice; nsv)
        U, Σ, V = B
    else
        U, Σ, V = svd(Data_ice)
    end

    # save in dictionary for later
    nfs = length(model_files)
    components_saved = joinpath(output_dir, "SVD_components_nfiles$(nfs)_r$(r).jld2")
    jldsave(components_saved; U, Σ, V, data_mean, data_ref)

    # prepare least square fit problem
    UΣ            = U*diagm(Σ)
    return UΣ, data_mean, data_ref, Σ, components_saved
end

function prepare_obs_SVD(gr, csv_dest, I_no_ocean, output_dir, fig_dir=""; input="dh_detrend")
    # retrieve standardized observations
    df_all = CSV.read(csv_dest, DataFrame)

    # use gdalgrid to rasterize atm point data
    df_atm      = df_all[df_all.source .== "atm", ["x", "y", input]]
    rename!(df_atm, input => "dh")    # just because gdalgrid is hard-wired here to look for "dh" variable
    tempname_csv = joinpath(output_dir, "dh_atm_temp_g$(gr).csv")
    tempname_nc  = joinpath(output_dir, "dh_atm_temp_g$(gr).nc")

    if !isfile(tempname_nc)
        CSV.write(tempname_csv, df_atm)
        gdalgrid(tempname_csv; gr, dest=tempname_nc)
    end
    obs_atm = NCDataset(tempname_nc)["Band1"][:]
    idx_atm = findall(.!ismissing.(obs_atm))

    # find aerodem indices
    id_df_aero = findall(df_all.source .== "aerodem")

    # merge everything in one matrix
    obs = zero(obs_atm)
    obs[df_all.idx[id_df_aero]] .= df_all[id_df_aero, input]
    obs[idx_atm]                .= obs_atm[idx_atm]

    # obtain I_obs and x_data vector
    I_obs      = findall(obs[I_no_ocean] .!= 0.0)
    x_data     = obs[I_no_ocean[I_obs]]

    if !isempty(fig_dir)
        # save matrix of observations as nc and plot
        save_netcdf(joinpath(dirname(fig_dir),"obs_all_gr$(gr).nc"), tempname_nc, [obs], ["dh"], Dict("dh"=>Dict{String,Any}()))
        Plots.heatmap(obs', clims=(-3,3), cmap=:bwr)
        Plots.savefig(joinpath(fig_dir, "obs_all.png"))
    end

    # remove temporary files
    # rm(tempname_csv)
    # rm(tempname_nc)

    return x_data, I_obs
end

function solve_optim(UΣ::Matrix{T}, I_obs::Vector{Int}, r::Int, λ::Real, x_data) where T <: Real
    @views A      = UΣ[I_obs,1:r]
    U_A, Σ_A, V_A = svd(A)
    D             = transpose(diagm(Σ_A))*diagm(Σ_A) + λ*I
    v_rec         = V_A * D^(-1) * transpose(diagm(Σ_A)) * transpose(U_A) * x_data
    x_rec         = UΣ[:,1:r]*v_rec
    return v_rec, x_rec
end

function SVD_reconstruction(λ::Real, r::Int, gr::Int, model_files::Vector{String}, csv_preproc, jld2_preproc, use_arpack=false)
    # define output paths
    main_output_dir = "output/SVD_reconstruction/"
    fig_dir = joinpath(main_output_dir, "figures")
    mkpath(main_output_dir)
    mkpath(fig_dir)

    # get reference DEM
    bedmachine_original, _ = create_bedmachine_grid(gr)
    _, ref_file            = create_grimpv2(gr, bedmachine_original)
    h_ref        = NCDataset(ref_file)["Band1"][:]
    h_ref        = replace_missing(h_ref, 0.0)

    # get I_no_ocean and (de-)standardization functions
    dict = load(jld2_preproc)
    @unpack I_no_ocean = dict
    standardize, destandardize = get_stddization_fcts(jld2_preproc)

    x_data, I_obs                          = prepare_obs_SVD(gr, csv_preproc, I_no_ocean, fig_dir)
    UΣ, data_mean, data_ref, _, saved_file = prepare_model(model_files, standardize, h_ref, I_no_ocean, r, main_output_dir; use_arpack) # read in model data and take svd to derive "eigen ice sheets"
    r                                      = min(size(UΣ,2), r)                                                                         # truncation of SVD cannot be higher than the second dimension of U*Σ
    v_rec, x_rec                           = solve_optim(UΣ, I_obs, r, λ, x_data)                                                       # derive analytical solution of regularized least squares

    # ad v_rec to output and save
    to_save = load(saved_file)
    to_save["v_rec"] = v_rec
    jldopen(saved_file, "w") do f
        [f[k] = m for (k,m) in to_save]
    end

    # calculate error and print
    nx, ny                  = size(NCDataset(ref_file)["Band1"])
    dif                     = zeros(F, nx,ny)
    dif[I_no_ocean[I_obs]] .= destandardize(x_rec[I_obs] .- x_data, h_ref[I_no_ocean[I_obs]], add_mean=false)
    err_mean                = mean(abs.(dif[I_no_ocean[I_obs]]))
    @printf("Mean absolute error: %1.1f m\n", err_mean)

    # retrieve matrix of reconstructed DEM
    dem_rec                 = zeros(F, nx,ny)
    dem_rec[I_no_ocean]     = data_ref .- destandardize(x_rec .+ data_mean, data_ref)
    dem_rec[dem_rec .== 0] .= no_data_value

    # save as nc file
    println("Saving file..")
    logλ        = Int(round(log(10, λ)))
    filename    = joinpath(main_output_dir,"rec_lambda_1e$logλ"*"_g$gr"*"_r$r.nc")
    layername   = "surface"
    attributes  = Dict(layername => Dict{String, Any}("long_name" => "ice surface elevation",
                                                      "standard_name" => "surface_altitude",
                                                      "units" => "m")
                        )
    save_netcdf(filename, ref_file, [dem_rec], [layername], attributes)

    # plot and save difference between reconstruction and observations
    save_netcdf(joinpath(main_output_dir, "dif_lambda_1e$logλ"*"_g$gr"*"_r$r.nc"), obs_aero_file, [dif], ["dif"], Dict("dif" => Dict{String,Any}()))
    Plots.heatmap(reshape(dif,nx,ny)', cmap=:bwr, clims=(-200,200), cbar_title="[m]", title="reconstructed - observations", size=(700,900))
    Plots.savefig(joinpath(fig_dir, "dif_lambda_1e$logλ"*"_g$gr"*"_r$r.png"))

    return filename, saved_file
end

function create_reconstructed_bedmachine(rec_file, bedmachine_file)
    # load datasets
    surfaceDEM = ncread(rec_file, "surface")
    bedDEM     = ncread(bedmachine_file, "bed")
    bedm_mask  = ncread(bedmachine_file, "mask")
    ice_mask   = (surfaceDEM .!= no_data_value) .&& (surfaceDEM .> bedDEM)

    # retrieve grid size
    x          = ncread(rec_file, "x")
    gr         = Int(x[2] - x[1])

    # calculate floating mask
    ρw            = 1030
    ρi            = 917
    Pw            = - ρw * bedDEM
    Pi            = ρi * (surfaceDEM - bedDEM)
    floating_mask = (Pw .> Pi) .&& (ice_mask)  # floating where water pressure > ice pressure at the bed

    # calculate mask
    new_mask = ones(eltype(bedm_mask), size(bedm_mask))
    new_mask[bedm_mask .== 0] .= 0 # ocean
    new_mask[ice_mask] .= 2
    new_mask[floating_mask]   .= 3

    # make sure DEM is zero everywhere on the ocean and equal to bed DEM on bedrock
    surfaceDEM[new_mask .== 0] .= no_data_value
    surfaceDEM[new_mask .== 1] .= bedDEM[new_mask .== 1]

    # calculate ice thickness
    h_ice = zeros(size(surfaceDEM)) .+ no_data_value
    h_ice[ice_mask] .= surfaceDEM[ice_mask] - bedDEM[ice_mask]
    h_ice[floating_mask] .= surfaceDEM[floating_mask] ./  (1-ρi/ρw)

    # save to netcdf file
    dest        = "output/bedmachine1980_reconstructed_g$(gr).nc"
    layers      = [surfaceDEM, bedDEM, h_ice, new_mask]
    layernames  = ["surface", "bed", "thickness", "mask"]
    template    = NCDataset(bedmachine_file)
    attributes  = get_attr(template, layernames)
    # overwrite some attributes
    sources_rec = Dict("surface"   => "svd reconstruction",
                       "bed"       => "Bedmachine-v5: Morlighem et al. (2022). IceBridge BedMachine Greenland, Version 5. Boulder, Colorado USA. NASA National Snow and Ice Data Center Distributed Active Archive Center. https://doi.org/10.5067/GMEVBWFLWA7X; projected on new grid with gdalwarp",
                       "thickness" => "computed from surface and bed",
                       "mask"      => "bedrock from Morlighem et al. (2022); ice, floating and ocean computed from surface and bed elevation"
                       )
    for l in layernames
        attributes[l]["source"] = sources_rec[l]
    end
    attributes["mask"]["long_name"] = "mask (0 = ocean, 1 = ice-free land, 2 = grounded ice, 3 = floating ice)"
    save_netcdf(dest, bedmachine_file, layers, layernames, attributes)

    return dest
end
