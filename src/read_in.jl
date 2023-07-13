using NetCDF, Glob, ProgressMeter

# function to read in model data
function read_model_data(;F::DataType=Float32,       # Float32 or Float64
                          which_files=nothing,       # indices of files used for training     ; e.g. 1:10, default all available
                          tsteps=nothing,            # indices of time steps used for training; e.g.  ""      ""
                          model_files)
    println("Reading in model data...")

    # determine indices for files and time steps
    if which_files === nothing
        which_files = 1:length(model_files)
    elseif which_files[end] > length(model_files)
        error("Number of files out of bound.")
    end
    files_out  = model_files[which_files]
    nf         = length(files_out)
    d          = ncread(files_out[1], "usurf")
    nx, ny, nt = size(d)
    if tsteps === nothing
        tsteps = 1:nt
    elseif tsteps[end] > nt
        error("Time steps out of bound.")
    end
    nt_out = length(tsteps)
    # build data matrix
    Data = zeros(F, nx*ny, nf*nt_out)
    @showprogress for (k, file) in enumerate(files_out)
        d = F.(ncread(file, "usurf")[:,:,tsteps])
        data = reshape(d, ny*nx, nt_out)
        Data[:,(k - 1 ) * nt_out + 1:k * nt_out] = data
    end
    return Data, nx, ny
end

"""
Get indices of cells with observations
"""
function get_indices(obs::Matrix{T}, mask_path::String, mask_name="Band1") where T<:Real
    # load imbie mask
    imbie_mask    = ncread(mask_path, mask_name)
    no_ocean_mask = findall((vec(obs) .> 0.0) .|| (vec(imbie_mask) .== 1))
    # get indices where there is data and ice, with respect to ice_mask
    R      = obs[no_ocean_mask]  # vector
    I_obs         = findall(R .> 0.0)
    return no_ocean_mask, I_obs
end
