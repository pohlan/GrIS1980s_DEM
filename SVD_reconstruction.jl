using svd_IceSheetDEM, NetCDF

# for running the script interactively
# ARGS = [
#         "--lambda", "1e7",
#         "--r", "300",
#         "--shp_file", "data/gris-imbie-1980/gris-outline-imbie-1980_updated.shp",
#         "--training_data", readdir("data/training_data/", join=true)...]

# retrieve command line arguments
parsed_args         = parse_commandline(ARGS)
training_data_files = parsed_args["training_data"]
outline_shp_file    = parsed_args["shp_file"]
use_arpack          = parsed_args["use_arpack"]
λ                   = parsed_args["λ"]
r                   = parsed_args["r"]

# make sure the training data set is not empty
@assert !isempty(training_data_files)
## choose a template file to make saving of netcdf files easier
template_file      = training_data_files[1]
## derive grid size in m from training data
x = ncread(template_file, "x")
const grd = Int(x[2] - x[1])   # assumes same grid size in both x and y direction

# make sure the shape file is available
if !isfile(outline_shp_file)
    error("shape file not found at " * outline_shp_file)
end

# Pre-process and standardize data #
csv_preprocessing, jld2_preprocessing, = prepare_obs(grd, outline_shp_file, nbins1=40, nbins2=50)

# reconstruction
rec_file, dict_file = SVD_reconstruction(λ, r, grd, training_data_files, csv_preprocessing, jld2_preprocessing; use_arpack)

# calculate the floating mask and create nc file according to the bedmachine template
dest        = joinpath("output", "reconstructions", "bedmachine1980_SVD_reconstruction_g$(Int(grd)).nc")
create_reconstructed_bedmachine(rec_file, dest)
