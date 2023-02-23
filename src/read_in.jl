# function to read in model data
function read_model_data(;which_files=nothing,       # indices of files used for training     ; e.g. 1:10, default all available
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
    Data = zeros(nx*ny, nf*nt_out)
    for (k, file) in enumerate(files_out)
        d = ncread(file, "usurf")[:,:,tsteps]
        data = reshape(d, ny*nx, nt_out)
        Data[:,(k - 1 ) * nt_out + 1:k * nt_out] = data
    end
    return Data, nx, ny
end
