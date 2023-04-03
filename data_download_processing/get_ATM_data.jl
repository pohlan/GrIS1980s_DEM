using DelimitedFiles, ProgressMeter, Glob, DataFrames, CSV, PyCall, Dates
import ArchGDAL as AG

path = "../data/"

## 1.) Download with python script from https://nsidc.org/data/data-access-tool/BLATM2/versions/1
pyinclude(fname) = (PyCall.pyeval_(read(fname, String), PyCall.pynamespace(Main), PyCall.pynamespace(Main), PyCall.Py_file_input, fname); nothing) # to be able to run an entire python script
pyinclude("nsidc-download_BLATM2.001_2023-03-19.py")
# move the files to a separate folder
mkpath(path*"ATM_raw")
fs = glob("BLATM2_*")
mv.(fs, path*"ATM_raw/".*fs, force=true)


## 2.) merge all flight files into one
files = glob(path*"ATM_raw/BLATM2_*_nadir2seg")
# appending DataFrames is much faster than doing e.g. 'd_final = [d_final; d_new]' or 'd_final = vcat(d_final, d_new)'
function merge_files(files)
    d_final = DataFrame()
    @showprogress for (n,f) in enumerate(files)
        startstring = "19"*split(f,"_")[3]
        starttime   = DateTime(startstring, dateformat"yyyymmdd")
        d           = readdlm(f)
        d_new       = DataFrame(d,:auto)
        d_new.x1    = starttime .+ Dates.Second.(Int.(round.(d[:,1])))
        d_new.x12   = d[:,3] .- 360
        if n == 1
            d_final = d_new
        else
            append!(d_final, d_new)
        end
    end
    return d_final
end
d_final = merge_files(files)
# save
CSV.write(path*"ATM_nadir2seg_all.csv", d_final)


## 3.) create regualar grid from scattered data using gdal_grid
### 3a) create a .vrt file with metadata for the csv file
open(path*"ATM_nadir2seg_all.vrt", "w") do io
    print(io,
"<OGRVRTDataSource>
    <OGRVRTLayer name=\"ATM_nadir2seg_all\">
        <SrcDataSource>"*path*"ATM_nadir2seg_all.csv</SrcDataSource>
        <SrcLayer>ATM_nadir2seg_all</SrcLayer>
        <LayerSRS>EPSG:4326</LayerSRS>
        <GeometryType>wkbPoint</GeometryType>
        <GeometryField encoding=\"PointFromColumns\" x=\"x12\" y=\"x2\" z=\"x4\"/>
    </OGRVRTLayer>
</OGRVRTDataSource>
")
end
### 3b) run gdal_grid
# nodata=-999.0 is important, otherwise gdalwarp will interpret the zeros as data
dataset   = AG.read(path*"ATM_nadir2seg_all.vrt")
options   = ["-a", "invdist:radius1=0.05:radius2=0.05:min_points=1:nodata=-9999.0",
             "-outsize", "1000", "1000",
             "-zfield", "x4",
             "-l", "ATM_nadir2seg_all"]
grid_data = AG.unsafe_gdalgrid(dataset, options)


## 4.) reproject to model grid using gdalwarp
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
for grid in [1200 1800]
    AG.unsafe_gdalwarp([grid_data], get_options(;grid); dest=path*"ATM_g$grid.nc")
end


## 5.) correct for geoid
##     bedmachine geoid = Eigen-6C4 Geoid - WGS84 Ellipsoid
##     ATM elevation    = elevation       - WGS84 Ellipsoid
##     ATM corrected    = ATM elevation   - bedmachine geoid

model_dataset = AG.read(path*"usurf_ex_gris_g1200m_v2023_RAGIS_id_0_1980-1-1_2020-1-1_YM.nc")
geoid_orig    = AG.read("NETCDF:"*path*"BedMachineGreenland-v5.nc:geoid")
for grid in [1200 1800]
    # reproject bedmachine geoid on model geometry
    geoid_file         = path*"bedm_geoid_g$grid.nc"
    AG.unsafe_gdalwarp([geoid_orig], get_options(;grid); dest=geoid_file)
    # apply correction
    geoid              = AG.read(AG.getband(AG.read(geoid_file),1))
    ATM                = AG.read(AG.getband(AG.read(path*"ATM_g$grid.nc"),1))
    ATM_corrected      = -9999.0 * ones(size(ATM))
    idx                = findall(ATM .!= 0)
    ATM_corrected[idx] = ATM[idx] - geoid[idx]

    # create netcdf file
    AG.create(
        path*"ATM_geoid_corrected_g$grid.nc",
        driver = AG.getdriver(model_dataset),
        width  = AG.width(model_dataset),
        height = AG.height(model_dataset),
        nbands = 1,
        dtype  = Float32
    ) do raster
        AG.write!(raster, ATM_corrected, 1)
        AG.setgeotransform!(raster, AG.getgeotransform(model_dataset))
        AG.setproj!(raster, AG.getproj(model_dataset))
    end
end
