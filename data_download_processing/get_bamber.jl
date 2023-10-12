# something seems a bit off with this DEM compared to GrIMP/model DEMs

using svd_IceSheetDEM, NetCDF

gr = 1200

# files
bamber_path = "data/bamber/"
bedm_file   = "data/bedmachine/bedmachine_g$(gr).nc"

# gdalwarp to bring it on right grid
bamb = svd_IceSheetDEM.gdalwarp(bamber_path*"GreenlandBamber1km_DEM.tif"; gr)  # or gdalwarp(path; gr, dest) to save it

# load geoid
geoid = ncread(bedm_file, "geoid")

# apply correction
bamber_corrected      = zeros(size(bamb))
idx                   = findall(bamb .!= 0)
bamber_corrected[idx] = bamb[idx] - geoid[idx]

# create netcdf file
dest = bamber_path*"bamber_geoid_corrected_g$gr.nc"
sample_path = "data/aerodem/aerodem_rm-filtered_geoid-corr_g$(gr).nc"
layername   = "surface"
attributes  = Dict(layername => Dict("long_name" => "ice surface elevation",
                                     "standard_name" => "surface_altitude",
                                     "units" => "m")
                   )
svd_IceSheetDEM.save_netcdf(dest, sample_path, [bamb[:,end:-1:1]], [layername], attributes)
