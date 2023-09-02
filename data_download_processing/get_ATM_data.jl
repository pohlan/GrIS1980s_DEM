using svd_IceSheetDEM, DelimitedFiles, ProgressMeter, Glob, DataFrames, CSV, PyCall, Dates
import ArchGDAL as AG

atm_path  = "data/ATM/"

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


## 3.) create regular grid from scattered data using gdal
# create a .vrt file with metadata for the csv file
open(atm_path*"ATM_nadir2seg_all.vrt", "w") do io
    print(io,
"<OGRVRTDataSource>
    <OGRVRTLayer name=\"ATM_nadir2seg_all\">
        <SrcDataSource>"*atm_path*"ATM_nadir2seg_all.csv</SrcDataSource>
        <SrcLayer>ATM_nadir2seg_all</SrcLayer>
        <LayerSRS>EPSG:4326</LayerSRS>
        <GeometryType>wkbPoint</GeometryType>
        <GeometryField encoding=\"PointFromColumns\" x=\"x12\" y=\"x2\" z=\"x4\"/>
    </OGRVRTLayer>
</OGRVRTDataSource>
")
end
# define options for gdalgrid
# nodata=-999.0 is important, otherwise gdalwarp will interpret the zeros as data
options   = ["-a", "invdist:radius1=0.05:radius2=0.05:min_points=1:nodata=-9999.0",
             "-outsize", "1000", "1000",
             "-zfield", "x4",
             "-l", "ATM_nadir2seg_all"]
# run gdalgrid and gdalwarp
gr      = 600
ATM = AG.read(atm_path*"ATM_nadir2seg_all.vrt") do source
    # from flight lines to grid
    AG.gdalgrid(source, options) do gridded
    # from arbitrary grid to model grid
        AG.gdalwarp([gridded],svd_IceSheetDEM.get_options(;gr); dest=atm_path*"ATM_g$gr.nc") do warped
            band = AG.getband(warped,1)
            AG.read(band)
        end
    end
end
# note about weird behaviour of ArchGDAL:
# the ATM matrix is upside down when plotted but will be flipped back when saving it to a netcdf further down


## 4.) correct for geoid
##     bedmachine geoid = Eigen-6C4 Geoid - WGS84 Ellipsoid
##     ATM elevation    = elevation       - WGS84 Ellipsoid
##     ATM corrected    = ATM elevation   - bedmachine geoid
bedm_path = "data/bedmachine/"
geoid              = ncread(bedm_path*"bedmachine_g$(gr).nc", "geoid")
ATM_corrected      = -9999.0 * ones(size(ATM))
idx                = findall(ATM .!= 0)
ATM_corrected[idx] = ATM[idx] - geoid[:,end:-1:1][idx]


## 5.) create netcdf file
training_path = "data/training_data_it0_$(gr)/"
sample_path   = training_path*"usurf_ex_gris_g$(gr)m_v2023_GIMP_id_0_1980-1-1_2020-1-1_YM.nc"
layername = "surface"
attributes  = Dict(layername => Dict("long_name" => "ice surface elevation",
                                     "standard_name" => "surface_altitude",
                                     "units" => "m")
                   )
svd_IceSheetDEM.save_netcdf(atm_path*"ATM_geoid_corrected_g$gr.nc", sample_path, [ATM_corrected[:,end:-1:1]], [layername], attributes)
