using svd_IceSheetDEM, NetCDF

# for running the script interactively
# ARGS = [
#         "--lambda", "1e5",
#         "--r", "377",
#         "--shp_file", "data/gris-imbie-1980/gris-outline-imbie-1980_updated.shp",
#         "--training_data", readdir("data/training_data_it0_1200", join=true)...]

parsed_args         = parse_commandline(ARGS)
training_data_files = parsed_args["training_data"]
outline_shp_file    = parsed_args["shp_file"]

# Make sure the training data set is not empty
@assert !isempty(training_data_files)
## choose a template file to make saving of netcdf files easier
template_file      = training_data_files[1]
## derive grid size in m from training data
x = ncread(template_file, "x")
const gr = Int(x[2] - x[1])   # assumes same grid size in both x and y direction

# Pre-process and standardize data #
csv_preprocessing, jld2_preprocessing, = prepare_obs(gr, outline_shp_file)


# -------------------------- #
# Part A: SVD reconstruction #
# -------------------------- #
# 1.) make sure the imbie shp file is available
if !isfile(outline_shp_file)
    error("shape file not found at " * outline_shp_file)
end

# 2.) run the svd solve_lsqfit

# retrieve command line arguments
λ           = parsed_args["λ"]     # regularization
r           = parsed_args["r"]
do_figures  = parsed_args["do_figures"]
use_arpack  = parsed_args["use_arpack"]
rec_file, dict_file = SVD_reconstruction(λ, r, gr, training_data_files, csv_preprocessing, jld2_preprocessing, use_arpack)

# 3.) calculate the floating mask and create nc file according to the bedmachine template
create_reconstructed_bedmachine(rec_file, bedmachine_file)  # ToDo --> after rf gneration??

# 4.) standardize residual, evaluate variogram and generate random fields
rf_files = SVD_random_fields(rec_file; n_fields=10)


# ------------------------------ #
# Part B: Interpolation approach #
# ------------------------------ #
grid_kriging = gr*2
interp_rec_file = geostats_interpolation(grid_kriging, gr, outline_shp_file, csv_preprocessing, jld2_preprocessing; maxn=10)
