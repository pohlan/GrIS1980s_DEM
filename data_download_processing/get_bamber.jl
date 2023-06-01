using svd_IceSheetDEM

grid = 1200

# gdalwarp to bring it on right grid
path = "data/GreenlandBamber1km_DEM.tif"
dest="data/bamber_test_g$grid.nc"
bamb = gdalwarp(path; grid)  # or gdalwarp(path; grid, dest) to save it

# load geoid
geoid = shortread("data/bedm_geoid_g$grid.nc")

# apply correction
bamber_corrected      = zeros(size(bamb))
idx                   = findall(bamb .!= 0)
bamber_corrected[idx] = bamb[idx] - geoid[idx]

# create netcdf file
dest = "data/bamber_geoid_corrected_g$grid.nc"
sample_path = "data/usurf_ex_gris_g1200m_v2023_RAGIS_id_0_1980-1-1_2020-1-1_YM.nc"
save_netcdf(bamber_corrected; dest, sample_path)
