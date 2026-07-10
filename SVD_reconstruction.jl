using GrIS1980s_DEM, NCDatasets, UnPack, JLD2

# for running the script interactively
# ARGS = [
#         "--lambda", "1e7",
#         "--r", "300",
#         "--shp_file", "data/outline/gris-outline-imbie-1980_updated.shp",
#         "--model_realizations", readdir("data/model_realizations/", join=true)...]

# retrieve command line arguments
parsed_args         = parse_commandline(ARGS)
model_realization_files = parsed_args["model_realizations"]
outline_shp_file    = parsed_args["shp_file"]
use_arpack          = parsed_args["use_arpack"]
λ                   = parsed_args["λ"]
r                   = parsed_args["r"]

# make sure the training data set is not empty
@assert !isempty(model_realization_files)
## choose a template file to make saving of netcdf files easier
template_file      = model_realization_files[1]
## derive grid size in m from training data
x = NCDataset(template_file)["x"][:]
const grd = Int(x[2] - x[1])   # assumes same grid size in both x and y direction

# make sure the shape file is available
if !isfile(outline_shp_file)
    error("shape file not found at " * outline_shp_file)
end

# Pre-process and standardize data #
csv_preprocessing, jld2_preprocessing, = prepare_obs(grd, outline_shp_file)

# reconstruction
rec_file, dict_file = SVD_reconstruction(λ, r, grd, model_realization_files, csv_preprocessing, jld2_preprocessing; use_arpack)

# calculate the floating mask and create nc file according to the bedmachine template
dest        = joinpath("output", "reconstructions", "rec_SVD_with_bedmachine_g$(Int(grd)).nc")
create_reconstructed_bedmachine(rec_file, dest)
