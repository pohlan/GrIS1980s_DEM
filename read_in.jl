# use historic image to find I (indices where there is data in 1980's)
# save as a JLD2 file
if !isfile("data/data1980s_indices.jld2")
    obs_file_historic ="data/aerodem_1978_1987_wgs84_g1800m.nc"
    obs_hist = ncread(obs_file_historic, "surface_altitude")
    nx, ny   = size(obs_hist)
    R_hist   = reshape(obs_hist, 1, length(obs_hist))
    I_marg   = findall(x->x!=-9999, R_hist[:])
    bool_intr   = obs_hist.==-9999
    for iy = 1:ny
        ix = 1
        while obs_hist[ix,iy] .== -9999
            bool_intr[ix,iy] = false
            ix += 1
            if ix > nx
                break
            end
        end
        ix = size(obs_hist,1)
        while obs_hist[ix,iy] .== -9999
            bool_intr[ix,iy] = false
            ix -= 1
            if ix < 1
                break
            end
        end
    end
    I_intr = findall(reshape(bool_intr,1,nx*ny)[:])
    jldsave("data/data1980s_indices.jld2"; I_marg, I_intr)
end

function read_model_data(;files::UnitRange{Int64},       # indices of files used for training
                          tsteps,      # indices of time steps used for training
                          model_files)
    println("Reading in model results...")
    training_files = model_files[files]
    nf = length(training_files)
    d = ncread(training_files[1], "usurf")
    nx, ny, nt = size(d)
    if length(tsteps) > nt error("Time steps out of bound.") end
    nttrain = length(tsteps)
    x = ncread(training_files[1], "x")
    y = ncread(training_files[1], "y")
    t = ncread(training_files[1], "time")
    # x-data: complete "data" sets of Greenland
    Data = zeros(nf * length(tsteps), nx*ny)
    for (k, training_file) in enumerate(training_files)
        d = ncread(training_file, "usurf")[:,:,tsteps]
        data = permutedims(reshape(d, ny*nx, nttrain))
        Data[(k - 1 ) * nttrain + 1:k * nttrain,:] = data
    end
    return Data, x, y, t, nx, ny
end
