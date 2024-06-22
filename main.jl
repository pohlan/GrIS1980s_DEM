using svd_IceSheetDEM, NetCDF

# for running the script interactively
# ARGS = [
#         "--lambda", "1e5",
#         "--r", "377",
#         "--shp_file", "data/gris-imbie-1980/gris-outline-imbie-1980_updated.shp",
#         "--training_data", readdir("data/training_data_it0_1200", join=true)...]

parsed_args         = parse_commandline(ARGS)
training_data_files = parsed_args["training_data"]
shp_file            = parsed_args["shp_file"]

# -------------------------- #
# Part A: SVD reconstruction #
# -------------------------- #

# 1.) make sure the training data set is not empty
@assert !isempty(training_data_files)
## choose a template file to make saving of netcdf files easier
template_file      = training_data_files[1]
## derive grid size in m from training data
x = ncread(template_file, "x")
const gr = Int(x[2] - x[1])   # assumes same grid size in both x and y direction

# 2.) move all the necessary bedmachine layers to the right grid (downloads the bedmachine-v5 if not available)
bedmachine_file = create_bedmachine_grid(gr, template_file)
bedmachine_path = splitdir(bedmachine_file)[1]

# 3.) make sure the imbie shp file is available
if !isfile(shp_file)
    error("shape file not found at " * shp_file)
end

# 4.) check if aerodem is available at the right grid, if not warp from available one or download/create from scratch
aerodem_g150, obs_file = create_aerodem(;gr, shp_file, bedmachine_path)

# 5.) get a netcdf mask from the imbie shp file
imbie_mask_file = create_imbie_mask(;gr, shp_file, sample_path=aerodem_g150)

# 6.) run the svd solve_lsqfit

# 377 -> findfirst(cumsum(Σ)./sum(Σ).>0.9)
# retrieve command line arguments
λ           = parsed_args["λ"]     # regularization
r           = parsed_args["r"]
do_figures  = parsed_args["do_figures"]
use_arpack  = parsed_args["use_arpack"]
rec_file    = SVD_reconstruction(λ, r, gr, imbie_mask_file, bedmachine_file, training_data_files[1:3], obs_file, do_figures, use_arpack)

# 5.) calculate the floating mask and create nc file according to the bedmachine template
create_reconstructed_bedmachine(rec_file, bedmachine_file)  # ToDo --> after rf gneration??


# ------------------------------------------------------------------------------- #
# Part B: residual analysis                                                       #
#   goal -> getting a distribution of reconstructed DEMs representing uncertainty #
# ------------------------------------------------------------------------------- #

# 1.) get ATM data
atm_file  = get_atm_file()
# 2.) get elevation change data from Sørensen et al., 2018
dh_obs_file, _   = create_dhdt_grid(;gr, startyr=1994, endyr=2010)
# 3.) standardize residual, evaluate variogram and generate random fields
rf_files = SVD_random_fields(rec_file, bedmachine_file, obs_file, atm_file, dh_obs_file, imbie_mask_file; n_fields=10)


# ------------------------------ #
# Part C: Interpolation approach #
# ------------------------------ #

interp_rec_file = geostats_interpolation(bedmachine_file, bedmachine_file, obs_file, atm_file, dh_obs_file, imbie_mask_file; nbins1=6, nbins2=12, maxn=50, n_fields=2, method=:sgs)
