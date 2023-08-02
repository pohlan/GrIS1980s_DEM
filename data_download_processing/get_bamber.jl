using svd_IceSheetDEM

gr = 1200

# gdalwarp to bring it on right grid
path = "data/GreenlandBamber1km_DEM.tif"
dest="data/bamber_test_g$gr.nc"
bamb = gdalwarp(path; gr)  # or gdalwarp(path; gr, dest) to save it

# load geoid
geoid = archgdal_read("data/bedm_geoid_g$gr.nc")

# apply correction
bamber_corrected      = zeros(size(bamb))
idx                   = findall(bamb .!= 0)
bamber_corrected[idx] = bamb[idx] - geoid[idx]

# create netcdf file
dest = "data/bamber_geoid_corrected_g$gr.nc"
sample_path = "data/usurf_ex_gris_g1200m_v2023_RAGIS_id_0_1980-1-1_2020-1-1_YM.nc"
save_netcdf(bamber_corrected; dest, sample_path)
