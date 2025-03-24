using svd_IceSheetDEM, NCDatasets, UnPack, JLD2

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
@unpack href_file = load(jld2_preprocessing)
@unpack standardize, destandardize = svd_IceSheetDEM.get_stddization_fcts(jld2_preprocessing)

# uncertainty estimation from cross-validation
cv_dict                    = load(get_cv_file_kriging(grd, maxn))
dem_ref                    = NCDataset(href_file)["surface"][:,:]
dif_destd                  = destandardize(cv_dict["difs"], cv_dict["binfield1"], cv_dict["h_ref"], add_mean=false)
dh_binned, bin_centers     = svd_IceSheetDEM.bin_equal_bin_size(cv_dict["h_ref"], dif_destd, 14)
sitp, std_uncertainty      = uncertainty_from_cv(dh_binned, bin_centers, dem_ref)
dest_file                  = get_std_uncrt_file(cv_dict["method"], grd)
svd_IceSheetDEM.save_netcdf(dest_file, href_file, [std_uncertainty], ["std_uncertainty"], Dict("std_uncertainty" => Dict{String,Any}()))

# reconstruction
tic = Base.time()
rec_file = geostats_interpolation(grd, outline_shp_file, csv_preprocessing, jld2_preprocessing; maxn)
toc = Base.time() - tic
dys  = round(toc ./ (3600*24), digits=1); println("Kriging took $dys days.")

# calculate the floating mask and create nc file according to the bedmachine template
dest        = joinpath("output", "reconstructions", "bedmachine1980_kriging_reconstruction_g$(Int(grd)).nc")
create_reconstructed_bedmachine(rec_file, dest)
