using svd_IceSheetDEM

# set target directories
main_output_dir = joinpath("output","validation")
mkpath(main_output_dir)
fig_dir         = joinpath(main_output_dir, "figures")
mkpath(fig_dir)

# for running the script interactively
# ARGS = [
#         "--shp_file", "data/gris-imbie-1980/gris-outline-imbie-1980_updated.shp",
#         "--grid_size", 600.0]

parsed_args         = parse_commandline(ARGS)
outline_shp_file    = parsed_args["shp_file"]
grd                 = parsed_args["grid_size"]
maxn                = parsed_args["maxn"]

# Pre-process and standardize data #
csv_preprocessing, jld2_preprocessing, = prepare_obs(grd, outline_shp_file, nbins1=40, nbins2=50)

# reconstruction
tic = Base.time()
rec_file = geostats_interpolation(grd, outline_shp_file, csv_preprocessing, jld2_preprocessing; maxn)
toc = Base.time() - tic
dys  = round(toc ./ (3600*24), digits=1); println("Kriging took $dys days.")

# calculate the floating mask and create nc file according to the bedmachine template
dest        = joinpath("output", "reconstructions", "bedmachine1980_kriging_reconstruction_g$(Int(grd)).nc")
create_reconstructed_bedmachine(rec_file, dest)
