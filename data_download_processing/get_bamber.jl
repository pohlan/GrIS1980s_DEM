import ArchGDAL as AG

bamber = AG.read("data/GreenlandBamber1km_DEM.tif")
bamb_data = AG.read(AG.getband(bamber, 1))

function get_options(;grid)
    x_min = - 678650
    y_min = -3371600
    x_max =   905350
    y_max = - 635600
    options = ["-overwrite",
               "-t_srs", "EPSG:3413",
               "-r", "average",
               "-co", "FORMAT=NC4",
               "-co", "COMPRESS=DEFLATE",
               "-co", "ZLEVEL=2",
               "-dstnodata", "0",
               "-te", "$x_min", "$y_min", "$x_max", "$y_max",
               "-tr", "$grid", "$grid"
               ]
    return options
end

grid = 1200
dest="data/bamber_g$grid.nc"
AG.unsafe_gdalwarp([bamber], get_options(;grid);dest)

geoid = AG.read(AG.getband(AG.read("data/bedm_geoid_g$grid.nc"),1))
bamb  = AG.read(AG.getband(AG.read(dest),1))

bamber_corrected      = zeros(size(bamb))
idx                   = findall(bamb .!= 0)
bamber_corrected[idx] = bamb[idx] - geoid[idx]

# create netcdf file
model_dataset = AG.read("data/usurf_ex_gris_g1200m_v2023_RAGIS_id_0_1980-1-1_2020-1-1_YM.nc")
AG.create(
    "data/bamber_geoid_corrected_g$grid.nc",
    driver = AG.getdriver(model_dataset),
    width  = AG.width(model_dataset),
    height = AG.height(model_dataset),
    nbands = 1,
    dtype  = Float32
) do raster
    AG.write!(raster, bamber_corrected, 1)
    AG.setgeotransform!(raster, AG.getgeotransform(model_dataset))
    AG.setproj!(raster, AG.getproj(model_dataset))
end
