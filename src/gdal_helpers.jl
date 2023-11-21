archgdal_read(file) = AG.read(AG.getband(AG.read(file),1))


"""
    get_options(; gr, cut_shp="", x_min=-678650, y_min=-3371600, x_max=905350, y_max=- 635600)
Returns a vector of Strings that can be used as an input for gdalwarp
"""
function get_options(;gr, cut_shp = "", srcnodata="", dstnodata="0", x_min = - 678650,
                                                                      y_min = -3371600,
                                                                      x_max =   905350,
                                                                      y_max = - 635600)
    options = ["-overwrite",
               "-t_srs", "EPSG:3413",
               "-r", "average",
               "-co", "FORMAT=NC4",
               "-co", "COMPRESS=DEFLATE",
               "-co", "ZLEVEL=2",
               "-te", "$x_min", "$y_min", "$x_max", "$y_max",  # y_max and y_min flipped otherwise GDAL is reversing the Y-dimension
               "-tr", "$gr", "$gr",
               ]
    if !isempty(cut_shp)
        append!(options, ["-cutline", cut_shp])
    end
    if !isempty(srcnodata)
        append!(options, ["-srcnodata", srcnodata])
    end
    if !isempty(dstnodata)
        append!(options, ["-dstnodata", dstnodata])
    end
    return options
end

"""
    gdalwarp(path; gr, kwargs...)

taken from https://discourse.julialang.org/t/help-request-using-archgdal-gdalwarp-for-resampling-a-raster-dataset/47039

# Example
```
out = gdalwarp("data/input.nc"; gr=1200, dest="data/output.nc")
```
"""
function gdalwarp(path::String; gr::Real, cut_shp="", srcnodata="", kwargs...)
    ds = AG.read(path) do source
        AG.gdalwarp([source], get_options(;gr, cut_shp, srcnodata); kwargs...) do warped
           band = AG.getband(warped, 1)
           AG.read(band)
       end
    end
    return ds
end
function gdalwarp(path::Vector{String}; gr::Real, cut_shp="", srcnodata="", kwargs...)
    source = AG.read.(path)
    ds = AG.gdalwarp(source, get_options(;gr, srcnodata)) do warped1
        if !isempty(cut_shp)
            AG.gdalwarp([warped1], get_options(;gr, cut_shp, srcnodata); kwargs...) do warped2
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

function get_attr(ds::NCDataset, layernames::Vector{String})
    # create dictionary of attributes
    attributes = Dict()
    for l in layernames
        var_attr = NCDatasets.OrderedDict()
        for att in filter(x -> !in(x, ["_FillValue"]), keys(ds[l].attrib))  # fillvalue can easily throw errors (https://docs.juliahub.com/NCDatasets/lxvtD/0.10.3/issues/#Defining-the-attributes-_FillValue,-add_offset,-scale_factor,-units-and-calendar)
            push!(var_attr, att => ds[l].attrib[att])
        end
        push!(attributes, l => var_attr)
    end
    return attributes
end

"""
    save_netcdf(A::Matrix; dest::String, sample_path::String)

- A:            Matrix of field to be saved
- dest:         destination where netcdf should be saved
- sample_path:  path of netcdf file that properties should be copied from
"""
function save_netcdf(dest::String, spatial_template_file::String, layers::Vector{T}, layernames::Vector{String}, attributes::Dict)  where T <: AbstractArray
    template  = NCDataset(spatial_template_file)
    crs_names = ["mapping", "polar_stereographic", "crs"]
    m    = crs_names[findfirst(in.(crs_names, (keys(template),)))]
    crs  = template[m]
    tx   = template["x"]
    ty   = template["y"]

    # make sure the template has the same size as layers
    nx, ny = size(layers[1])
    @assert (nx, ny) == (length(tx), length(ty))

    # write fields and their attributes
    rm(dest, force=true) # NCDataset(..., "c") cannot be done twice in the same Julia session (Permission denied error); if attempted, the file gets corrupted and needs to be deleted
    ds  = NCDataset(dest, "c")
    defDim(ds, "x", nx)
    defDim(ds, "y", ny)
    for (field, data) in zip(layernames, layers)
        if !haskey(attributes[field], "grid_mapping")       # without the grid_mapping attribute gdalwarp doesn't work!
            push!(attributes[field], "grid_mapping" => m)
        end
        defVar(ds, field, data, ("x", "y"), attrib = attributes[field])
    end
    defVar(ds, m, Char, (), attrib = crs.attrib)
    defVar(ds, "x", tx[:], ("x",), attrib = tx.attrib)
    defVar(ds, "y", ty[:], ("y",), attrib = ty.attrib)

    # global attributes
    gr = template["x"][2] - template["x"][1]
    ds.attrib["spacing in m"] = "$gr"

    close(ds)
end

"""
    Do gdalwarp of mask manually because I couldn't find an option in gdalwarp to make the search radius larger
    Necessary because bedrock mask needs to be conservatively large not leaving too many individual ice cells in between.
"""
function mask_downsamle(bedmachine_fullres, templ, r)
    bedm_src = ncread(bedmachine_fullres, "mask")
    x_src    = ncread(bedmachine_fullres, "x")
    y_src    = ncread(bedmachine_fullres, "y")
    x_dst    = ncread(templ, "x")
    y_dst    = ncread(templ, "y")
    bedm_dst = zeros(length(x_dst), length(y_dst))
    # nx, ny = size(bedm_src)
    @showprogress for iy in axes(bedm_dst, 2)
        for ix in axes(bedm_dst, 1)
            ix_src = findall(abs.(x_dst[ix] .- x_src) .< r)
            iy_src = findall(abs.(y_dst[iy] .- y_src) .< r)
            if !any(isempty.([ix_src, iy_src]))
                # println("$ix, $iy")
                pts_in_range = bedm_src[ix_src, iy_src]
                if sum(pts_in_range .== 0) < 0.5*length(pts_in_range)
                    pts_in_range[pts_in_range .== 0] .= 10
                end
                bedm_dst[ix,iy] = minimum(pts_in_range)
            else
                bedm_dst[ix,iy] = 0
            end
        end
    end
    return bedm_dst
end


function create_bedmachine_grid(gr, spatial_template_file, rad_mask_downsmpl=1000)
    bedmachine_path = "data/bedmachine/"
    dest_file = bedmachine_path * "bedmachine_g$(gr).nc"
    # if file exists already, do nothing
    if isfile(dest_file)
        return dest_file
    end

    # check if the original bedmachine is there already, otherwise download
    mkpath(bedmachine_path)
    bedmachine_original = bedmachine_path*"BedMachineGreenland-v5.nc"
    if !isfile(bedmachine_original)
        println("Downloading bedmachine..")
        url_bedm = "https://n5eil01u.ecs.nsidc.org/ICEBRIDGE/IDBMG4.005/1993.01.01/BedMachineGreenland-v5.nc"
        Downloads.download(url_bedm, bedmachine_original)
    end

    # gdalwarp
    println("Using gdalwarp to project bedmachine on model grid..")
    layernames = ["mask", "geoid", "bed", "surface", "thickness"]
    fieldvals  = Matrix[]
    for fn in layernames
        if fn == "mask"
            new = mask_downsamle(bedmachine_original, spatial_template_file, rad_mask_downsmpl)   # 1000 is the radius over which the minimum value is taken
        else
            new = gdalwarp("NETCDF:"*bedmachine_original*":"*fn;  gr)[:,end:-1:1]
        end
        push!(fieldvals, new)  # flip y-axis due to behaviour of ArchGDAL
    end

    # save as netcdf
    attr_template = NCDataset(bedmachine_original)
    attributes = get_attr(attr_template, layernames)
    save_netcdf(dest_file, spatial_template_file, fieldvals, layernames, attributes)
    return dest_file
end

"""
copied from https://gist.github.com/scls19fr/9ea2fd021d5dd9a97271da317bff6533
"""
function get_table_from_html(input::AbstractString)
    r = read(Downloads.download(input), String)      # recognizes .netcr file automatically, in contrast to HTTP request
    h = Gumbo.parsehtml(r)
    @assert strip(nodeText(eachmatch(Selector("title"), h.root)[1])) != "Earthdata Login" "Make sure .netrc file is available and correct."
    qs = eachmatch(Selector("table"), h.root)
    tables = []
    for helm_table in qs
        column_names = String[]
        d_table = OrderedDict{String, Vector{String}}()
        for (i, row) in enumerate(eachmatch(Selector("tr"), helm_table))
            if (i == 1)
                for (j, colh) in enumerate(eachmatch(Selector("th"), row))
                    colh_text = strip(nodeText(colh))
                    while (colh_text in column_names)  # column header must be unique
                        colh_text = colh_text * "_2"
                    end
                    push!(column_names, colh_text)
                end
            else
                if (i == 2)
                    for colname in column_names
                        d_table[colname] = Vector{String}()
                    end
                end
                for (j, col) in enumerate(eachmatch(Selector("td"), row))
                    col_text = strip(nodeText(col))
                    colname = column_names[j]
                    push!(d_table[colname], col_text)
                end
            end
        end
        df = DataFrame(d_table)
        if length(qs) == 1
            return df
        end
        push!(tables, df)
    end
    return tables
end

function create_aerodem(;gr, shp_file, bedmachine_path, kw="")
    aerodem_path = "data/aerodem/"
    get_aero_file(gr) = aerodem_path * "aerodem_rm-filtered_geoid-corr_g$(gr).nc"
    aerodem_g150_file = get_aero_file(150)
    aerodem_gr_file   = get_aero_file(gr)

    if isfile(aerodem_gr_file)
        return aerodem_g150_file, aerodem_gr_file
    elseif !isfile(aerodem_g150_file)

        # create aerodem, for some reason the cutting with the shapefile outline only works for smaller grids
        # otherwise GDALError (CE_Failure, code 1): Cutline polygon is invalid.
        raw_path  = aerodem_path*"raw/"
        mkpath(raw_path)

        # download
        if isempty(readdir(raw_path))
            println("Downloading aerodem tif files, this may take a few minutes...")
            url_DEMs          = "https://www.nodc.noaa.gov/archive/arc0088/0145405/1.1/data/0-data/G150AERODEM/DEM/"
            url_reliablt_mask = "https://www.nodc.noaa.gov/archive/arc0088/0145405/1.1/data/0-data/G150AERODEM/ReliabilityMask/"
            for url in [url_DEMs, url_reliablt_mask]
                df = get_table_from_html(url)
                tif_files = df.Name[endswith.(df.Name, ".tif") .&& .!occursin.("carey", df.Name) .&& occursin.(kw, df.Name)]
                missing_files = tif_files[.!isfile.(raw_path.*tif_files)]
                Downloads.download.(url .* missing_files, raw_path .* missing_files)
            end
        end

        aerodem_files = glob(raw_path * "aerodem_*.tif")
        rm_files      = glob(raw_path * "rm_*.tif")

        # gdalwarp
        println("Using gdalwarp to merge the aerodem mosaics into one DEM, taking a while..")
        merged_aero_dest = aerodem_path*"merged_aerodem_g150.nc"
        aero_rm_filt     = gdalwarp(aerodem_files; gr=150, cut_shp=shp_file, dest=merged_aero_dest)
        rel_mask         = gdalwarp(     rm_files; gr=150, cut_shp=shp_file)

        # filter for observations where reliability value is low
        aero_rm_filt[rel_mask .< 40] .= 0.0  # only keep values with reliability of at least xx

        # apply geoid correction
        geoid = gdalwarp("NETCDF:"*joinpath(bedmachine_path,"BedMachineGreenland-v5.nc:geoid");  gr=150)
        idx                                = findall(aero_rm_filt .!= 0.0)
        aero_rm_filt[idx]                .-= geoid[idx]
        aero_rm_filt[aero_rm_filt .< 0.0] .= 0.0

        # save as netcdf
        sample_path = merged_aero_dest
        layername   = "surface"
        attributes  = Dict(layername => Dict("long_name" => "ice surface elevation",
                                             "standard_name" => "surface_altitude",
                                             "units" => "m")
                            )
        save_netcdf(aerodem_g150_file, sample_path, [aero_rm_filt[:,end:-1:1]], [layername], attributes)
    end

    # gdalwarp to desired grid
    gdalwarp(aerodem_g150_file; gr, srcnodata="0.0", dest=aerodem_gr_file)

    return aerodem_g150_file, aerodem_gr_file
end

function create_imbie_mask(;gr, shp_file, sample_path)
    imbie_path = "data/gris-imbie-1980/"
    imbie_mask_file = imbie_path * "imbie_mask_g$(gr).nc"
    # if file exists already, do nothing
    if isfile(imbie_mask_file)
        return imbie_mask_file
    end

    # a bit ugly and cumbersome, due to the unability to resolve two issues:
    # 1) cut_shp doesn't work wih too large grids and if src file doesn't correspond to grid size
    # 2) the ArchGDAL gdalwarp function doesn't take Matrices as an input (or at least I haven't figured out how)
    #    need to give a filename as argument
    println("Creating imbie mask netcdf..")
    mkpath(imbie_path)
    sample = archgdal_read(sample_path)
    ones_m = ones(size(sample))
    fname_ones = "temp1.nc"
    fname_mask = "temp2.nc"
    layername   = "mask"
    attributes = Dict(layername => Dict())
    save_netcdf(fname_ones, sample_path, [ones_m], [layername], attributes)
    gdalwarp(fname_ones; gr=150, cut_shp=shp_file, dest=fname_mask)
    gdalwarp(fname_mask; gr, dest=imbie_mask_file)
    rm(fname_ones, force=true)
    rm(fname_mask, force=true)
    return imbie_mask_file
end

function get_nc_from_flightlines(pt_data::DataFrame, bedm_file::String, dest_file::String; correct_geoid=false)
    # load coordinates of grid to project on
    x      = ncread(bedm_file, "x")
    y      = ncread(bedm_file, "y")
    dx     = x[2]-x[1]
    dy     = y[2]-y[1]
    nx, ny = length(x), length(y)

    # transform atm coordinates
    coords = [[pt_data.lon[i], pt_data.lat[i]] for i in eachindex(pt_data.lat)]
    coords_proj = AG.reproject(coords, ProjString("+proj=longlat +datum=WGS84 +no_defs"), EPSG(3413))  # hard-coded source and target coordinate systems
    x_proj = first.(coords_proj)
    y_proj = last.(coords_proj)

    # reproject on grid using averages weighted by inverse distance
    println("Projecting flight line values on model grid...")
    r_inv     = zeros(Float32, nx, ny)
    grid_data = zeros(Float32, nx, ny)
    npts      = zeros(Int, nx, ny)
    n = 1
    for (xp, yp, zp) in zip(x_proj, y_proj, pt_data.z)
        rx, ix = findmin(abs.(xp .- x))
        ry, iy = findmin(abs.(yp .- y))
        r = sqrt(rx^2 + ry^2)
        if abs(xp - x[ix]) < 0.5dx && abs(yp - y[iy]) < 0.5dy
            npts[ix,iy] += 1
            r_inv[ix,iy] += 1 / r
            grid_data[ix,iy] += zp / r
        end
        n += 1
    end
    i_z            = grid_data .!= 0
    grid_data[i_z] = grid_data[i_z] ./ r_inv[i_z]

    # correct for geoid
    if correct_geoid
        geoid = ncread(bedm_file, "geoid")
        grid_data[i_z] -= geoid[i_z]
    end

    # save as netcdf file
    svd_IceSheetDEM.save_netcdf(dest_file, bedm_file, [grid_data], ["surface"], Dict("surface" => Dict()))
    return
end

function create_atm_grid(gr, bedm_file::String, kw="")
    atm_path = "data/ATM/"
    atm_file = atm_path*"ATM_elevation_geoid_corrected_g$(gr).nc"
    # if file exists already, do nothing
    if isfile(atm_file)
        return atm_file
    end

    # download files
    raw_path  = atm_path*"ATM_raw/"
    mkpath(raw_path)
    if isempty(readdir(raw_path))
        println("Donloading ATM elevation data...")
        atm_url = "https://n5eil01u.ecs.nsidc.org/PRE_OIB/BLATM2.001/"
        tb1 = get_table_from_html(atm_url)
        folders = tb1.Name[(startswith.(tb1.Name, "1993") .|| startswith.(tb1.Name, "1994")) .&& occursin.(kw, tb1.Name)]
        for fld in folders
            df = get_table_from_html(atm_url*fld)
            files_to_download = df.Name[endswith.(df.Name, "_nadir2seg")]
            missing_files = files_to_download[.!isfile.(raw_path.*files_to_download)]    # only download missing files in case some are downloaded already
            Downloads.download.(atm_url*fld.*missing_files, raw_path.*missing_files)
        end
    end
    files     = glob(raw_path*"BLATM2_*_nadir2seg")

    # merge all flight files into one
    # appending DataFrames is much faster than doing e.g. 'd_final = [d_final; d_new]' or 'd_final = vcat(d_final, d_new)'
    function merge_files(files)
        d_final = DataFrame()
        @showprogress for (n,f) in enumerate(files)
            startstring = "19"*split(basename(f),"_")[2]
            starttime   = DateTime(startstring, dateformat"yyyymmdd")
            d           = readdlm(f)
            time        = starttime .+ Dates.Second.(Int.(round.(d[:,1])))
            lon         = d[:,3]
            lat         = d[:,2]
            z           = d[:,4]
            d_new       = DataFrame((;time, lon, lat, z))
            if n == 1
                d_final = d_new
            else
                append!(d_final, d_new)
            end
        end
        return d_final
    end
    df = merge_files(files)

    get_nc_from_flightlines(df, bedm_file, atm_file; correct_geoid=true)
    return atm_file
end
