import ArchGDAL as AG

shortread(file) = AG.read(AG.getband(AG.read(file),1))

"""
    get_options(; grid, cut_shp="", x_min=-678650, y_min=-3371600, x_max=905350, y_max=- 635600)
Returns a vector of Strings that can be used as an input for gdalwarp
"""
function get_options(;grid, cut_shp = "", srcnodata="", x_min = - 678650,
                                                        y_min = -3371600,
                                                        x_max =   905350,
                                                        y_max = - 635600)
    options = ["-overwrite",
               "-t_srs", "EPSG:3413",
               "-r", "average",
               "-co", "FORMAT=NC4",
               "-co", "COMPRESS=DEFLATE",
               "-co", "ZLEVEL=2",
               "-dstnodata", "0",
               "-te", "$x_min", "$y_min", "$x_max", "$y_max",  # y_max and y_min flipped otherwise GDAL is reversing the Y-dimension
               "-tr", "$grid", "$grid",
               ]
    if !isempty(cut_shp)
        append!(options, ["-cutline", cut_shp])
    end
    if !isempty(srcnodata)
        append!(options, ["-srcnodata", srcnodata])
    end
    return options
end

"""
    gdalwarp(path; grid, kwargs...)

taken from https://discourse.julialang.org/t/help-request-using-archgdal-gdalwarp-for-resampling-a-raster-dataset/47039

# Example
```
out = gdalwarp("data/input.nc"; grid=1200, dest="data/output.nc")
```
"""
function gdalwarp(path::String; grid::Real, cut_shp="", srcnodata="", kwargs...)
    ds = AG.read(path) do source
        AG.gdalwarp([source], get_options(;grid, cut_shp, srcnodata); kwargs...) do warped
           band = AG.getband(warped, 1)
           AG.read(band)
       end
    end
    return ds
end
function gdalwarp(path::Vector{String}; grid::Real, cut_shp="", srcnodata="", kwargs...)
    source = AG.read.(path)
    ds = AG.gdalwarp(source, get_options(;grid, srcnodata)) do warped1
        if !isempty(cut_shp)
            AG.gdalwarp([warped1], get_options(;grid, cut_shp, srcnodata); kwargs...) do warped2
                band = AG.getband(warped2, 1)
                AG.read(band)
            end
        else
            band = AG.getband(warped1, 1)
            AG.read(band)
        end
    end
    return ds
end

"""
    save_netcdf(A::Matrix; dest::String, sample_path::String)

- A:            Matrix of field to be saved
- dest:         destination where netcdf should be saved
- sample_path:  path of netcdf file that properties should be copied from
"""
function save_netcdf(A::Matrix; dest::String, sample_path::String, )
    sample_dataset = AG.read(sample_path)
    AG.create(
        dest,
        driver = AG.getdriver(sample_dataset),
        width  = AG.width(sample_dataset),
        height = AG.height(sample_dataset),
        nbands = 1,
        dtype  = Float32,
        options = ["-co", "COMPRESS=DEFLATE"] # reduces the file size
    ) do raster
        AG.write!(raster, A, 1)
        AG.setgeotransform!(raster, AG.getgeotransform(sample_dataset))
        AG.setproj!(raster, AG.getproj(sample_dataset))
    end
end
