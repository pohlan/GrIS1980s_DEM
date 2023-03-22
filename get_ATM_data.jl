using DelimitedFiles, ProgressMeter, Glob, DataFrames, CSV, PyCall
import ArchGDAL as AG

## 1.) Download with python script from https://nsidc.org/data/data-access-tool/BLATM2/versions/1
pyinclude(fname) = (PyCall.pyeval_(read(fname, String), PyCall.pynamespace(Main), PyCall.pynamespace(Main), PyCall.Py_file_input, fname); nothing) # to be able to run an entire python script
pyinclude("nsidc-download_BLATM2.001_2023-03-19.py")
# move the files to a separate folder
mkpath("data/ATM")
fs = glob("BLATM2_*")
mv.(fs, "data/ATM/".*fs, force=true)


## 2.) merge all flight files into one
files = glob("data/ATM/BLATM2_*_nadir2seg")
# appending DataFrames is much faster than doing e.g. 'd_final = [d_final; d_new]' or 'd_final = vcat(d_final, d_new)'
function merge_files(files)
    d_final = DataFrame()
    @showprogress for (n,f) in enumerate(files)
        startstring      = "19"*split(f,"_")[2]
        starttime        = DateTime(startstring, dateformat"yyyymmdd")
        d                = readdlm(f)
        d_new            = DataFrame(d,:auto)
        d_new.x1         = starttime .+ Dates.Second.(round.(d[:,1]))
        d_new.alt_xcoord = d[:,3] .- 360
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
CSV.write("data/ATM_nadir2seg_all.csv", d_final)


## 3.) create regualar grid from scattered data using gdal_grid
## nodata=-999.0 is important, otherwise gdalwarp will interpret the zeros as data
dataset   = AG.read("data/ATM_nadir2seg_all.vrt")
options   = ["-a", "invdist:radius1=0.05:radius2=0.05:min_points=1:nodata=-999.0",
             "-outsize", "111", "222",
             "-zfield", "x4",
             "-l", "ATM_nadir2seg_all"]
grid_data = AG.unsafe_gdalgrid(dataset, options)
# band1 = AG.getband(grid_data, 1)
# b_matrix = AG.read(band1)


## 4.) reproject to model grid using gdalwarp
x_min=-678650
y_min=-3371600
x_max=905350
y_max=-635600
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

for grid in [1200 1800]
    AG.unsafe_gdalwarp([a], options; dest="data/ATM_g$grid.nc")
end
