#!/bin/bash
# chmod u+x DEM_alignment_GDAL.sh
set -x -e
export HDF5_USE_FILE_LOCKING=FALSE
options='-overwrite  -t_srs EPSG:3413 -r average -co FORMAT=NC4 -co COMPRESS=DEFLATE -co ZLEVEL=2 '
x_min=-678650
y_min=-3371600
x_max=905350
y_max=-635600
data=data/bambgrl_dem_5km_corrected.tif
# data=data/aerodem_g1800m_geoid_corrected_1978_1987_mean.nc
# data=NETCDF:data/BedMachineGreenland-v5.nc:surface
for grid in 1800; do
    datpism=data/bamber_gris_g${grid}.nc
    gdalwarp $options  -dstnodata 0 -te $x_min $y_min $x_max $y_max -tr $grid $grid  $data $datpism
    ncrename -v Band1,surface_altitude datpism
done
