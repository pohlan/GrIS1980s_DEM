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

function prepare_model(model_files, I_no_ocean, r, dest_file=""; use_arpack=false)
    # load model data and calculate difference to reference DEM
    Data_ice  = read_model_data(;model_files,I_no_ocean)

    # compute SVD
    if use_arpack
        nsv = min(r, size(Data_ice, 2)-1) # actual truncation is later, but takes too long if r is unnecessarily high here
        B, _ = svds(Data_ice; nsv)
        U, Σ, V = B
    else
        U, Σ, V = svd(Data_ice)
    end

    # save in dictionary for later
    if !isempty(dest_file)
        nfiles = length(model_files)
        jldsave(dest_file; U, Σ, V, nfiles)
    end

    # prepare least square fit problem
    UΣ            = U*diagm(Σ)
    return UΣ, Σ
end

function prepare_obs_SVD(grd, csv_dest, I_no_ocean, output_dir)
    # retrieve standardized observations
    df_all = CSV.read(csv_dest, DataFrame)

    # use gdalgrid to rasterize atm point data
    df_atm      = df_all[df_all.source .== "atm", ["x", "y", "h"]]
    rename!(df_atm, "h" => "dh")    # just because gdalgrid is hard-wired here to look for "dh" variable
    tempname_csv = joinpath(output_dir, "atm_temp_g$(grd).csv")
    tempname_nc  = joinpath(output_dir, "atm_temp_g$(grd).nc")

    if !isfile(tempname_nc)
        CSV.write(tempname_csv, df_atm)
        gdalgrid(tempname_csv; grd, dest=tempname_nc)
    end
    obs_atm = NCDataset(tempname_nc)["Band1"][:,:]
    idx_atm = findall(.!ismissing.(obs_atm))

    # find aerodem indices
    id_df_aero = findall(df_all.source .== "aerodem")
    @assert !isempty(id_df_aero)

    # merge everything in one matrix
    obs = zero(obs_atm)
    obs[df_all.idx[id_df_aero]] .= df_all[id_df_aero, "h"]
    obs[idx_atm]                .= obs_atm[idx_atm]

    # obtain I_obs and x_data vector
    I_obs      = findall(obs[I_no_ocean] .!= 0.0)
    x_data     = obs[I_no_ocean[I_obs]]

    # save matrix of observations as netcdf
    obs[obs .== 0] .= no_data_value
    save_netcdf(joinpath(output_dir,"obs_all_gr$(grd).nc"), tempname_nc, [obs], ["h"], Dict("h"=>Dict{String,Any}()))
    return x_data, I_obs
end

function solve_optim(UΣ::Matrix{T}, I_obs::Vector{Int}, r::Int, λ::Real, x_data) where T <: Real
    @views A      = UΣ[I_obs,1:r]
    v_rec         = (transpose(A)*A + λ*I)^(-1) * transpose(A) * x_data
    # v_rec         = (transpose(A)*W*A + λ*I)^(-1) * transpose(A)*W*x_data
    x_rec         = UΣ[:,1:r]*v_rec
    return v_rec, x_rec
end

function SVD_reconstruction(λ::Real, r::Int, grd::Int, model_files::Vector{String}, csv_preproc, jld2_preproc; use_arpack=false)
    # define output paths
    main_output_dir = joinpath("output","reconstructions")
    mkpath(main_output_dir)

    # get I_no_ocean and (de-)standardization functions
    dict = load(jld2_preproc)
    @unpack I_no_ocean, href_file = dict

    h_ref        = NCDataset(href_file)["surface"]
    # rel_mask     = nomissing(NCDataset("data/aerodem/rm_g$(grd).nc")["Band1"][:,:], 1)

    saved_file    = joinpath(main_output_dir, "SVD_components_g$(grd)_rec.jld2")
    UΣ, _         = prepare_model(model_files, I_no_ocean, r, saved_file; use_arpack) # read in model data and take svd to derive "eigen ice sheets"
    x_data, I_obs = prepare_obs_SVD(grd, csv_preproc, I_no_ocean, main_output_dir)
    r             = min(size(UΣ,2), r)                                                                         # truncation of SVD cannot be higher than the second dimension of U*Σ
    # W = Diagonal(rel_mask[I_no_ocean[I_obs]])
    v_rec, x_rec                           = solve_optim(UΣ, I_obs, r, λ, x_data)                                                       # derive analytical solution of regularized least squares

    # ad v_rec to output and save
    to_save = load(saved_file)
    to_save["v_rec"] = v_rec
    jldopen(saved_file, "w") do f
        [f[k] = m for (k,m) in to_save]
    end

    # calculate error and print
    nx, ny                  = size(h_ref)
    dif                     = zeros(F, nx,ny)
    dif[I_no_ocean[I_obs]] .= x_rec[I_obs] .- x_data
    err_mean                = mean(abs.(dif[I_no_ocean[I_obs]]))
    @printf("Mean absolute error: %1.1f m\n", err_mean)

    # retrieve matrix of reconstructed DEM
    dem_rec                 = zeros(F, nx,ny)
    dem_rec[I_no_ocean]    .= x_rec
    dem_rec[dem_rec .<= 0] .= no_data_value

    # estimated error from cross-validation
    std_uncertainty = NCDataset(get_std_uncrt_file("SVD", grd))["std_uncertainty"][:,:]

    # save as nc file
    println("Saving file..")
    logλ        = Int(round(log(10, λ)))
    filename    = joinpath(main_output_dir,"rec_g$(grd)_lambda_1e$(logλ)_r$(r).nc")
    attributes  = Dict("surface" => Dict{String, Any}("long_name" => "ice surface elevation",
                                                      "units" => "m"),
                       "std_uncertainty" => Dict{String,Any}("long_name" => "standard deviation of error estimated from cross-validation",
                                                                 "units" => "m")
                        )
    save_netcdf(filename, href_file, [dem_rec, std_uncertainty], ["surface", "std_uncertainty"], attributes)

    # save difference between reconstruction and observations
    save_netcdf(joinpath(main_output_dir, "dif_lambda_1e$logλ"*"_g$grd"*"_r$r.nc"), href_file, [dif], ["dif"], Dict("dif" => Dict{String,Any}()))

    return filename, saved_file
end
