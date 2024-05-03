# function to read in model data
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
    return Data, nx, ny
end

"""
Get indices of cells with observations
"""
function get_indices(obs::Matrix{T}, imbie_path::String, bedm_path::String, mask_name="Band1") where T<:Real
    # load imbie mask
    imbie_mask = ncread(imbie_path, mask_name)
    grimp_mask = ncread(bedm_path, "mask")
    ice_mask   = findall( (vec(grimp_mask) .!= 1) .&& vec(imbie_mask) .> 0.0)
    # get indices where there is data and ice, with respect to ice_mask
    R          = obs[ice_mask]  # vector
    I_obs      = findall(R .> 0.0)
    return ice_mask, I_obs
end
