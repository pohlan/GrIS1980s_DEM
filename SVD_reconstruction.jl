using svd_IceSheetDEM, NCDatasets, UnPack, JLD2

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
x = NCDataset(template_file)["x"][:]
const grd = Int(x[2] - x[1])   # assumes same grid size in both x and y direction

# make sure the shape file is available
if !isfile(outline_shp_file)
    error("shape file not found at " * outline_shp_file)
end

# Pre-process and standardize data #
csv_preprocessing, jld2_preprocessing, = prepare_obs(grd, outline_shp_file, nbins1=40, nbins2=50)
@unpack href_file = load(jld2_preprocessing)

# uncertainty estimation from cross-validation
cv_dict                = load(get_cv_file_SVD(grd, length(training_data_files)))
ir, iλ                 = findfirst(cv_dict["rs"].==r), findfirst(cv_dict["λs"].==λ)
dem_ref                = NCDataset(href_file)["surface"][:,:]
dh_binned, bin_centers = svd_IceSheetDEM.bin_equal_bin_size(dem_ref[cv_dict["idx"]], cv_dict["m_difs"][iλ,ir], 14)
sitp, std_uncertainty  = uncertainty_from_cv(dh_binned, bin_centers, dem_ref)
dest_file              = get_std_uncrt_file(cv_dict["method"], grd)
svd_IceSheetDEM.save_netcdf(dest_file, href_file, [std_uncertainty], ["std_uncertainty"], Dict("std_uncertainty" => Dict{String,Any}()))

# reconstruction
rec_file, dict_file = SVD_reconstruction(λ, r, grd, training_data_files, csv_preprocessing, jld2_preprocessing; use_arpack)

# calculate the floating mask and create nc file according to the bedmachine template
dest        = joinpath("output", "reconstructions", "bedmachine1980_SVD_reconstruction_g$(Int(grd)).nc")
create_reconstructed_bedmachine(rec_file, dest)
