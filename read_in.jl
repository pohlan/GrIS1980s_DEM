# 40 time steps in total, 9 files


function make_training_data(;rtrain::UnitRange{Int64},     # indices of files used for training
                            tsteps::UnitRange{Int64},      # indices of time steps used for training
                            t_y::Int64                     # timestep used for fitting y-data in I
)
    println("Reading in data...")
    model_files = glob("data/usurf_gris_g1800m_v5_RAGIS_id_*_1980-1-1_2020-1-1_YM.nc")
    #fend = length(model_files)÷(1/rtrain)
    training_files = model_files[rtrain]
    nf = length(training_files)
    d = ncread(training_files[1], "usurf")
    nx, ny, nt = size(d)
    if length(tsteps) > nt error("Time steps out of bound.") end
    nttrain = length(tsteps)
    x = ncread(training_files[1], "x")
    y = ncread(training_files[1], "y")

    # x-data: complete "data" sets of Greenland
    Data = zeros(nf * length(tsteps), nx*ny)
    for (k, training_file) in enumerate(training_files)
        d = ncread(training_file, "usurf")[:,:,tsteps]
        data = transpose(reshape(d, ny * nx, nttrain))
        Data[(k - 1 ) * nttrain + 1:k * nttrain, :] = data
    end

    # use historic image to find I (indices where there is data in 1980's)
    obs_file_historic ="data/aerodem_1978_1987_wgs84_g1800m.nc"
    obs_hist = ncread(obs_file_historic, "surface_altitude")
    R_hist = reshape(obs_hist, 1, nx * ny)
    I = findall(x->x!=-9999, R_hist)
    I_inv = findall(x->x==-9999, R_hist)

    # y-data:
    obs_artif = ncread(model_files[end],"usurf")[:,:,t_y]
    R = reshape(obs_artif, 1, nx * ny);

    return Data, I, I_inv, R
end