# Assumes that the name of the files are all fixed

using svd_IceSheetDEM, NetCDF

const F = Float32 # Julia default is Float64 but that kills the process for the full training data set if r is too large

# for running the script interactively
# ARGS = [
#         "--lambda", "1e5",
#         "--r", "377",
#         "--imbie_shp_file", "data/gris-imbie-1980/gris-outline-imbie-1980_updated.shp",
#         "--training_data", readdir("data/training_data_it0_1200", join=true)...]

parsed_args = parse_commandline(ARGS)

imbie_path          = "data/gris-imbie-1980/"
aerodem_path        = "data/aerodem/"
bedmachine_path     = "data/bedmachine/"
training_data_files = parsed_args["training_data"]
imbie_shp_file      = parsed_args["imbie_shp_file"]

# 1.) make sure the training data set is not empty
@assert !isempty(training_data_files)
## choose a template file to make saving of netcdf files easier
template_file      = training_data_files[1]
## derive grid size in m from training data
x = ncread(template_file, "x")
const gr = Int(x[2] - x[1])   # assumes same grid size in both x and y direction

# 2.) move all the necessary bedmachine layers to the right grid (downloads the bedmachine-v5 if not available)
bedmachine_file = bedmachine_path * "bedmachine_g$(gr).nc"
if !isfile(bedmachine_file)
    create_bedmachine_grid(gr, bedmachine_path, template_file)
end

# 3.) make sure the imbie shp file is available
if isnothing(imbie_shp_file)
    error("no imbie shape file provided")
elseif !isfile(imbie_shp_file)
    error("imbie shape file not found at " * imbie_shp_file)
end

# 4.) check if aerodem is available at the right grid, if not warp from available one or download/create from scratch
aerodem_g150 = aerodem_path * "aerodem_rm-filtered_geoid-corr_g150.nc"
obs_file     = aerodem_path*"aerodem_rm-filtered_geoid-corr_g$(gr).nc"
if !isfile(aerodem_path * "aerodem_rm-filtered_geoid-corr_g$(gr).nc")
    if !isfile(aerodem_g150)
        # create aerodem, for some reason the cutting with the shapefile outline only works for smaller grids
        # otherwise GDALError (CE_Failure, code 1): Cutline polygon is invalid.
        create_aerodem(aerodem_path, imbie_shp_file, bedmachine_path)
    end
    gdalwarp(aerodem_g150; gr, srcnodata="0.0", dest=obs_file)
end

# 5.) get a netcdf mask from the imbie shp file
imbie_mask = imbie_path * "imbie_mask_g$(gr).nc"
if !isfile(imbie_mask)
    create_imbie_mask(gr; imbie_path, imbie_shp_file, sample_path=aerodem_g150)
end


# 6.) run the svd solve_lsqfit

# 377 -> findfirst(cumsum(Σ)./sum(Σ).>0.9)
# retrieve command line arguments
λ           = F(parsed_args["λ"])     # regularization
r           = parsed_args["r"]
do_figure   = parsed_args["do_figure"]
rec_file    = solve_lsqfit(F, λ, r, gr, imbie_mask, training_data_files, obs_file, do_figure)

# 5.) calculate the floating mask
create_reconstructed_bedmachine(rec_file, bedmachine_file)
