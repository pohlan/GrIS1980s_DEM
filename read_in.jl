# use historic image to find I (indices where there is data in 1980's)
# save as a JLD2 file
if !isfile("data/data1980s_indices.jld2")
    obs_file_historic ="data/aerodem_1978_1987_wgs84_g1800m.nc"
    obs_hist = ncread(obs_file_historic, "surface_altitude")
    R_hist = reshape(obs_hist, length(obs_hist), 1)
    I = findall(x->x!=-9999, R_hist)
    I_inv = findall(x->x==-9999, R_hist)
    jldsave("data/data1980s_indices.jld2"; I, I_inv)
end

function xdata_training(;files::UnitRange{Int64},       # indices of files used for training
                         tsteps::UnitRange{Int64},      # indices of time steps used for training
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
    # x-data: complete "data" sets of Greenland
    Data = zeros(nx*ny,nf * length(tsteps))
    for (k, training_file) in enumerate(training_files)
        d = ncread(training_file, "usurf")[:,:,tsteps]
        data = reshape(d, ny*nx, nttrain)
        Data[:,(k - 1 ) * nttrain + 1:k * nttrain] = data
    end
    return Data, nx, ny
end
