using GrIS1980s_DEM, UnPack, JLD2, NCDatasets

parsed_args         = parse_commandline(ARGS)
outline_shp_file    = parsed_args["shp_file"]
grd                 = parsed_args["grid_size"]

# Pre-processing
csv_preprocessing, jld2_preprocessing = GrIS1980s_DEM.prepare_obs(grd, outline_shp_file)

# interpolation
tic = Base.time()
rec_file = geostats_interpolation(grd, outline_shp_file, csv_preprocessing, jld2_preprocessing, ℓ_block=5e4, δl=1.2e5)
toc = Base.time() - tic
dys = round(toc ./ (3600*24), digits=1); println("GP took $dys days.")

# calculate the floating mask and create nc file according to the bedmachine template
dest = joinpath("output", "reconstructions", "bedmachine1980_GP_reconstruction_g$(Int(grd)).nc")
create_reconstructed_bedmachine(rec_file, dest)
