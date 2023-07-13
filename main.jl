# Assumes that the name of the files are all fixed

using svd_IceSheetDEM

imbie_path         = "data/gris-imbie-1980/"
aerodem_path       = "data/aerodem/"
bedmachine_path    = "data/bedmachine/"
training_data_path = "data/training_data_it0_1200/"

# 1.) make sure the training data set is not empty
@assert !isempty(training_data_path)
## choose a template file to make saving of netcdf files easier
template_path      = training_data_path*readdir(training_data_path)[1]
## derive grid size in m from training data
const gr = parse(Int, split(template_path, "_")[end-6][2:end-1])

# 2.) move all the necessary bedmachine layers to the right grid (downloads the bedmachine-v5 if not available)
bedmachine_file = bedmachine_path * "bedmachine_g$(gr).nc"
if !isfile(bedmachine_file)
    create_bedmachine_grid(gr; bedmachine_path, template_path)
end

# 3.) check if aerodem is available at the right grid, if not warp from available one or download/create from scratch
aerodem_g150 = aerodem_path * "aerodem_rm-filtered_geoid-corr_g150.nc"
obs_file     = aerodem_path*"aerodem_rm-filtered_geoid-corr_g$(gr).nc"
if !isfile(aerodem_path * "aerodem_rm-filtered_geoid-corr_g$(gr).nc")
    if !isfile(aerodem_g150)
        # create aerodem, for some reason the cutting with the shapefile outline only works for smaller grids
        # otherwise GDALError (CE_Failure, code 1): Cutline polygon is invalid.
        create_aerodem(aerodem_path, imbie_path, bedmachine_path)
    end
    gdalwarp(aerodem_g150; grid=gr, srcnodata="0.0", dest=obs_file)
end

# 4.) make sure that the imbie shp file is downloaded and get a netcdf mask of the right grid
imbie_mask = imbie_path * "imbie_mask_g$(gr).nc"
if !isfile(imbie_path * "gris-outline-imbie-1980.shp")
    @error "shape file of imbie outline not downloaded"
end
if !isfile(imbie_mask)
    create_imbie_mask(gr; imbie_path, sample_path=aerodem_g150)
end


# 4.) run the svd solve_lsqfit
ARGS = ["--save",
        "--lambda", "1e5",
        "--r", "377"]
# 377 -> findfirst(cumsum(Σ)./sum(Σ).>0.9)
include("SVD_lsqfit.jl")
rec_file = solve_it(ARGS, gr, imbie_mask, training_data_path, obs_file)

# 5.) calculate the floating mask
include("floating_mask.jl")
create_reconstructed_bedmachine(obs_file, rec_file, bedmachine_file)
