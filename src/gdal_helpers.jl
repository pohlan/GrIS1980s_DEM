archgdal_read(file) = AG.read(AG.getband(AG.read(file),1))

"""
Defines the domain containing the Greenland ice sheet that all grids are projected on, including the projection
"""
@kwdef struct DomainLimits{T}
    x_min::T = - 678650
    y_min::T = -3371600
    x_max::T =   905350
    y_max::T = - 635600
    crs::String = "EPSG:3413"
end

"""
    get_options(; gr, cut_shp="", x_min=-678650, y_min=-3371600, x_max=905350, y_max=- 635600)
Returns a vector of Strings that can be used as an input for gdalwarp
"""
function gdalwarp_options(;gr, cut_shp = "", srcnodata="", resampling="bilinear")
    lims = DomainLimits()
    options = ["-overwrite",
               "-t_srs", lims.crs,
               "-r", resampling,
               "-co", "FORMAT=NC4",
               "-co", "COMPRESS=DEFLATE",
               "-co", "ZLEVEL=2",
               "-dstnodata", string(no_data_value),
               "-te", string(lims.x_min), string(lims.y_min), string(lims.x_max), string(lims.y_max),
               "-tr", "$gr", "$gr",
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
    gdalwarp(path; gr, kwargs...)

taken from https://discourse.julialang.org/t/help-request-using-archgdal-gdalwarp-for-resampling-a-raster-dataset/47039

# Example
```
out = gdalwarp("data/input.nc"; gr=1200, dest="data/output.nc")
```
"""
function gdalwarp(path::String; gr::Real, cut_shp="", srcnodata="", resampling="bilinear", kwargs...)
    ds = AG.read(path) do source
        AG.gdalwarp([source], gdalwarp_options(;gr, cut_shp, srcnodata, resampling); kwargs...) do warped
           band = AG.getband(warped, 1)
           AG.read(band)
       end
    end
    return ds
end
function gdalwarp(path::Vector{String}; gr::Real, cut_shp="", srcnodata="", resampling="bilinear", kwargs...)
    source = AG.read.(path)
    ds = AG.gdalwarp(source, gdalwarp_options(;gr, srcnodata, resampling); kwargs...) do warped1
        if !isempty(cut_shp)
            AG.gdalwarp([warped1], gdalwarp_options(;gr, cut_shp, srcnodata, resampling); kwargs...) do warped2
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

function gdalgrid(filenm::String; gr::Real, kwargs...)
    filenm_noext  = splitext(filenm)[1]             # including path but without extension
    filenm_basenm = splitext(basename(filenm))[1]   # no path, no extension

    # get domain limits and crs of target grid
    lims = DomainLimits()

    # create .vrt file containing metadata for csv file
    open(filenm_noext*".vrt", "w") do io
        print(io,
    "<OGRVRTDataSource>
        <OGRVRTLayer name=\""*filenm_basenm*"\">
            <SrcDataSource>"*filenm*"</SrcDataSource>
            <SrcLayer>"*filenm_basenm*"</SrcLayer>
            <LayerSRS>"*lims.crs*"</LayerSRS>
            <GeometryType>wkbPoint</GeometryType>
            <GeometryField encoding=\"PointFromColumns\" x=\"x\" y=\"y\" z=\"dh\"/>
        </OGRVRTLayer>
    </OGRVRTDataSource>
    ")
    end

    # define options for gdalgrid
    rad     = 2*gr
    options = ["-a", "invdist:radius=$rad:min_points=1:nodata=$(string(no_data_value))",
               "-co", "FORMAT=NC4",
               "-co", "COMPRESS=DEFLATE",
               "-co", "ZLEVEL=2",
               "-txe", string(lims.x_min), string(lims.x_max),
               "-tye", string(lims.y_min), string(lims.y_max),
               "-tr", "$gr", "$gr"]
    display(options)
    # call gdalgrid
    ds = AG.read(filenm_noext*".vrt") do source
        AG.gdalgrid(source, options; kwargs...) do gridded
            band = AG.getband(gridded, 1)
            AG.read(band)
        end
    end
    return
end

function get_attr(ds::NCDataset, layernames::Vector{String})
    # create dictionary of attributes
    attributes = Dict{String, Dict{String, Any}}()
    for l in layernames
        var_attr = NCDatasets.OrderedDict()
        for att in filter(x -> !in(x, ["_FillValue"]), keys(ds[l].attrib))  # fillvalue is prescribed, not necessarily the same as in template
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
function save_netcdf(dest::String, spatial_template_file::String, layers::Vector{T}, layernames::Vector{String}, attributes::Dict{String, Dict{String, Any}})  where T <: AbstractArray
    template  = NCDataset(spatial_template_file)
    crs_names = ["mapping", "polar_stereographic", "crs", "spatial_ref"]
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
        # make sure the grid mapping attribute is there and has the correct entry
        if !haskey(attributes[field], "grid_mapping")       # without the grid_mapping attribute gdalwarp doesn't work!
            push!(attributes[field], "grid_mapping" => m)
        else
            attributes[field]["grid_mapping"] = m
        end
        if field == "mask"
            defVar(ds, field, data, ("x", "y"), attrib = attributes[field])
        else
            defVar(ds, field, data, ("x", "y"), attrib = attributes[field], fillvalue = no_data_value)
        end
    end
    defVar(ds, m, Char, (), attrib = crs.attrib)
    defVar(ds, "x", tx[:], ("x",), attrib = tx.attrib)
    defVar(ds, "y", ty[:], ("y",), attrib = ty.attrib)

    # global attributes
    gr = template["x"][2] - template["x"][1]
    ds.attrib["spacing in m"] = "$gr"

    close(ds)
end

# """
#     Do gdalwarp of mask manually because I couldn't find an option in gdalwarp to make the search radius larger
#     Necessary because bedrock mask needs to be conservatively large not leaving too many individual ice cells in between.
# """
# function mask_downsample(bedmachine_fullres, templ, r)
#     bedm_src = ncread(bedmachine_fullres, "mask")
#     x_src    = ncread(bedmachine_fullres, "x")
#     y_src    = ncread(bedmachine_fullres, "y")
#     x_dst    = ncread(templ, "x")
#     y_dst    = ncread(templ, "y")
#     bedm_dst = zeros(Int8, length(x_dst), length(y_dst))
#     @showprogress for iy in axes(bedm_dst, 2)
#         for ix in axes(bedm_dst, 1)
#             ix_src = findall(abs.(x_dst[ix] .- x_src) .< r)
#             iy_src = findall(abs.(y_dst[iy] .- y_src) .< r)
#             if !any(isempty.([ix_src, iy_src]))
#                 pts_in_range = bedm_src[ix_src, iy_src]
#                 if sum(pts_in_range .== 0) < 0.5*length(pts_in_range)
#                     pts_in_range[pts_in_range .== 0] .= 10
#                 end
#                 bedm_dst[ix,iy] = minimum(pts_in_range)
#             end
#         end
#     end
#     return bedm_dst
# end


function create_bedmachine_grid(gr)
    bedmachine_path = "data/bedmachine/"
    dest_file = joinpath(bedmachine_path, "bedmachine_g$(gr).nc")
    bedmachine_original = bedmachine_path*"BedMachineGreenland-v5.nc"
    # if file exists already, do nothing
    if isfile(dest_file)
        return bedmachine_original, dest_file
    end

    # check if the original bedmachine is there already, otherwise download
    mkpath(bedmachine_path)
    if !isfile(bedmachine_original)
        println("Downloading bedmachine..")
        url_bedm = "https://n5eil01u.ecs.nsidc.org/ICEBRIDGE/IDBMG4.005/1993.01.01/BedMachineGreenland-v5.nc"
        Downloads.download(url_bedm, bedmachine_original)
    end

    # gdalwarp
    println("Using gdalwarp to project bedmachine on model grid..")
    layernames = ["geoid", "bed", "surface", "thickness", "mask"]
    fieldvals  = Matrix[]
    attr_template = NCDataset(bedmachine_original)
    spatial_template_file = "template_file.nc"
    for fn in layernames
        resampling = fn == "mask" ? "max" : "bilinear"
        new = gdalwarp("NETCDF:"*bedmachine_original*":"*fn;  gr, srcnodata=string(attr_template.attrib["no_data"]), resampling, dest=spatial_template_file)[:,end:-1:1]
        push!(fieldvals, new)  # flip y-axis due to behaviour of ArchGDAL
    end

    # save as netcdf
    attributes = get_attr(attr_template, layernames)
    save_netcdf(dest_file, spatial_template_file, fieldvals, layernames, attributes)
    rm(spatial_template_file)
    return bedmachine_original, dest_file
end

function create_grimpv2(gr, bedmachine_original; kw="")
    data_path = joinpath("data","grimpv2")
    get_dest_file(gr) = joinpath(data_path, "grimpv2_geoid_corrected_g$(gr).nc")
    dest_g150_file    = get_dest_file(150)
    dest_gr_file      = get_dest_file(gr)

    if isfile(dest_gr_file)
        return dest_g150_file, dest_gr_file
    elseif !isfile(dest_g150_file)

        # set up paths
        raw_path = joinpath(data_path, "raw")
        mkpath(data_path)
        mkpath(raw_path)

        # download
        if isempty(readdir(raw_path))
            println("Download grimp v2 surface DEM...")
            url_grimpv2 = "https://n5eil01u.ecs.nsidc.org/MEASURES/NSIDC-0715.002/2008.05.15/"
            t = get_table_from_html(url_grimpv2)
            tif_files = t.Name[endswith.(t.Name, ".tif") .&& occursin.("_dem_", t.Name) .&& [any(contains.(name, kw)) for name in t.Name]]
            Downloads.download.(joinpath.(url_grimpv2, tif_files), joinpath.(raw_path,tif_files))
        end

        # gdalwarp
        grimp_files = readdir(raw_path, join=true)
        filter!(endswith(".tif"),grimp_files)
        merged_g150 = joinpath(data_path,"merged_grimpv2_g150.nc")
        println("Using gdalwarp to merge the grimpv2 mosaics into one DEM, taking a while..")
        grimp_g150  = gdalwarp(grimp_files; gr=150, srcnodata="0.0", dest=merged_g150)
        grimp_g150  = grimp_g150[:,end:-1:1]

        # apply geoid correction
        geoid = gdalwarp("NETCDF:"*bedmachine_original*":geoid";  gr=150)[:,end:-1:1]
        idx   = findall(grimp_g150 .!= no_data_value .&& .!ismissing.(geoid))
        grimp_g150[idx] .-= geoid[idx]
        grimp_g150[grimp_g150 .< 0.0] .= no_data_value

        # save as netcdf
        sample_path = merged_g150
        save_netcdf(dest_g150_file, sample_path, [grimp_g150], ["surface"], Dict("surface"=>Dict{String,Any}()))
        rm(merged_g150)
    end

    # gdalwarp to desired grid
    gdalwarp(dest_g150_file; gr, srcnodata=string(no_data_value), dest=dest_gr_file)

    return dest_g150_file, dest_gr_file
end

function get_surface_file(ref1, bedm_file; remove_geoid=false)
    ds_ref1 = NCDataset(ref1)
    surf    = ds_ref1["Band1"][:]
    if remove_geoid  # when comparing to elevation data that is not geoid corrected
        geoid    = NCDataset(bedm_file)["geoid"][:]
        surf   .+= geoid
    end
    surf[ismissing.(surf)] .= no_data_value
    verticalrefname = remove_geoid ? "ellipsoid" : "geoid"
    out_name = splitext(ref1)[1]*"_surface_"*verticalrefname*splitext(ref1)[2]
    save_netcdf(out_name, ref1, [Float32.(surf)], ["surface"], Dict("surface" => Dict{String,Any}()))
    return out_name
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

function create_aerodem(gr, outline_shp_file, bedmachine_original, reference_file_g150; kw="")
    aerodem_path      = joinpath("data","aerodem")
    get_aero_file(gr) = joinpath(aerodem_path, "aerodem_rm-filtered_geoid-corr_g$(gr).nc")
    # get_rm_file(gr)   = aerodem_path * "rm_g$(gr).nc"
    aerodem_g150_file = get_aero_file(150)
    aerodem_gr_file   = get_aero_file(gr)
    # rm_g150_file      = get_rm_file(150)
    # rm_gr_file        = get_rm_file(gr)

    if isfile(aerodem_gr_file) #&& isfile(rm_gr_file)
        return aerodem_g150_file, aerodem_gr_file #, rm_gr_file
    elseif !isfile(aerodem_g150_file)

        # create aerodem, for some reason the cutting with the shapefile outline only works for smaller grids
        # otherwise GDALError (CE_Failure, code 1): Cutline polygon is invalid.
        raw_path  = joinpath(aerodem_path,"raw")
        fig_path  = joinpath(aerodem_path, "figures")
        mkpath(aerodem_path)
        mkpath(raw_path)
        mkpath(fig_path)

        # download
        if isempty(readdir(raw_path))
            println("Downloading aerodem tif files, this may take a few minutes...")
            url_DEMs          = "https://www.nodc.noaa.gov/archive/arc0088/0145405/1.1/data/0-data/G150AERODEM/DEM/"
            url_reliablt_mask = "https://www.nodc.noaa.gov/archive/arc0088/0145405/1.1/data/0-data/G150AERODEM/ReliabilityMask/"
            for url in [url_DEMs, url_reliablt_mask]
                df = get_table_from_html(url)
                tif_files = df.Name[endswith.(df.Name, ".tif") .&& .!occursin.("carey", df.Name) .&& occursin.(kw, df.Name)]
                missing_files = tif_files[.!isfile.(raw_path.*tif_files)]
                Downloads.download.(url .* missing_files, joinpath.(raw_path, missing_files))
            end
        end

        aerodem_files = glob(joinpath(raw_path, "aerodem_*.tif"))
        rm_files      = glob(joinpath(raw_path, "rm_*.tif"))
        # gdalwarp
        merged_aero_dest = joinpath(aerodem_path, "merged_aerodem_g150.nc")
        merged_rm_dest   = joinpath(aerodem_path, "merged_rm_g150.nc")
        println("Using gdalwarp to merge the aerodem mosaics into one DEM, taking a while..")
        aero_rm_filt     = gdalwarp(aerodem_files; gr=150, srcnodata=string(0.0), dest=merged_aero_dest)
        rel_mask         = gdalwarp(     rm_files; gr=150, srcnodata=string(0.0), dest=merged_rm_dest)

        # filter for observations where reliability value is low
        aero_rm_filt[rel_mask .< 40] .= no_data_value  # only keep values with reliability of at least xx

        # apply geoid correction
        geoid = gdalwarp("NETCDF:"*bedmachine_original*":geoid";  gr=150)
        idx                                = findall(aero_rm_filt .!= no_data_value)
        aero_rm_filt[idx]                .-= geoid[idx]
        aero_rm_filt[aero_rm_filt .< 0.0] .= no_data_value

        # save as netcdf
        sample_path = merged_aero_dest
        layername   = "surface"
        attributes  = Dict(layername => Dict{String, Any}("long_name"     => "ice surface elevation",
                                                          "standard_name" => "surface_altitude",
                                                          "units"         => "m")
                            )
        not_aligned_file = splitext(aerodem_g150_file)[1]*"_not_aligned"*splitext(aerodem_g150_file)[2]
        save_netcdf(not_aligned_file, sample_path, [aero_rm_filt[:,end:-1:1]], [layername], attributes)
        # save_netcdf(rm_g150_file, sample_path, [Float32.(rel_mask[:,end:-1:1])], ["reliability mask"], Dict("reliability mask" => Dict{String, Any}()))

        # co-registration
        fig_name = joinpath(fig_path, "coregistration_before_after.jpg")
        py_coregistration(reference_file_g150, not_aligned_file, aerodem_g150_file, outline_shp_file, fig_name) # note: output is upside down but will be reversed in gdalwarp below

        rm(merged_aero_dest)
        rm(merged_rm_dest)
    end

    # gdalwarp to desired grid
    gdalwarp(aerodem_g150_file; gr, cut_shp=outline_shp_file, srcnodata=string(no_data_value), dest=aerodem_gr_file)
    # gdalwarp(rm_g150_file; gr, srcnodata=string(no_data_value), dest=rm_gr_file)

    return aerodem_g150_file, aerodem_gr_file
end

function create_outline_mask(gr, outline_shp_file, sample_file)
    outline_path = "data/gris-imbie-1980/"
    outline_mask_file = joinpath(outline_path, "outline_mask_g$(gr).nc")
    # if file exists already, do nothing
    if isfile(outline_mask_file)
        return outline_mask_file
    end

    # a bit ugly because
    # the ArchGDAL gdalwarp function doesn't take matrices as an input (or at least I haven't figured out how)
    # --> need to save an nc file with ones and give the filename as an input
    println("Creating outline mask netcdf..")
    mkpath(outline_path)
    sample = archgdal_read(sample_file)
    ones_m = ones(size(sample))
    fname_ones = "temp1.nc"
    layername   = "mask"
    attributes = Dict(layername => Dict{String, Any}())
    save_netcdf(fname_ones, sample_file, [ones_m], [layername], attributes)
    gdalwarp(fname_ones; gr, cut_shp=outline_shp_file, dest=outline_mask_file)
    rm(fname_ones, force=true)
    return outline_mask_file
end

function get_atm_raw_file(kw)
    atm_path = "data/ATM/"
    atm_file = joinpath(atm_path,"ATM_nadir2seg_all.csv")
    # if file exists already, do nothing
    if isfile(atm_file)
        return atm_file
    end

    # download files
    raw_path  = joinpath(atm_path,"raw/")
    mkpath(raw_path)
    if isempty(readdir(raw_path))
        println("Downloading ATM elevation data...")
        atm_url = "https://n5eil01u.ecs.nsidc.org/PRE_OIB/BLATM2.001/"
        tb1 = get_table_from_html(atm_url)
        folders = tb1.Name[(startswith.(tb1.Name, "1993") .|| startswith.(tb1.Name, "1994")) .&& [any(contains.(name, kw)) for name in tb1.Name]]
        for fld in folders
            df = get_table_from_html(atm_url*fld)
            files_to_download = df.Name[endswith.(df.Name, "_nadir2seg")]
            missing_files = files_to_download[.!isfile.(raw_path.*files_to_download)]    # only download missing files in case some are downloaded already
            Downloads.download.(atm_url*fld.*missing_files, joinpath.(raw_path,missing_files))
        end
    end
    files     = glob(joinpath(raw_path,"BLATM2_*_nadir2seg"))

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
    CSV.write(atm_file, df)
    return atm_file
end

function get_atm_dh_file(ref_file, bedm_file, blockspacing; kw="")
    raw_atm_data_file = get_atm_raw_file(kw)

    atm_dh_dest_file = joinpath(dirname(raw_atm_data_file), "dh_ref_minus_atm_blocksize_$(blockspacing)m.csv")
    # if file exists already, do nothing
    if isfile(atm_dh_dest_file)
        return atm_dh_dest_file
    end

    # need a file with only one band for interp_points operation in xdem python;
    # plus need to reference to ellipsoid to be consistent with atm before differencing the two
    ref_surface_ellipsoid = get_surface_file(ref_file, bedm_file, remove_geoid=true)
    ref_surface_geoid     = get_surface_file(ref_file, bedm_file, remove_geoid=false)
    # interpolate GrIMP DEM on ATM points and calculate difference
    dh_interpolated_file = joinpath(dirname(raw_atm_data_file), "atm_dh_interpolated.csv")
    py_point_interp(ref_surface_ellipsoid, ref_surface_geoid, raw_atm_data_file, dh_interpolated_file)

    # block reduce with python package verde; average data that is heavily oversampled in direction of flight
    py_block_reduce(dh_interpolated_file, atm_dh_dest_file, blockspacing)

    # delete all temporary files
    rm(ref_surface_ellipsoid)
    rm(ref_surface_geoid)
    rm(dh_interpolated_file)

    return atm_dh_dest_file
end

function get_atm_grid(gr, ref_file, bedm_file, blockspacing=gr)
    atm_dh_file   = get_atm_dh_file(ref_file, bedm_file, blockspacing)
    println(atm_dh_file)
    atm_grid_file = joinpath(dirname(atm_dh_file), "dh_ref_minus_atm_g$(gr).nc")
    if isfile(atm_grid_file)
        return atm_grid_file
    end
    println("Use gdalgrid to project atm flightline data on grid...")
    gdalgrid(atm_dh_file; gr, dest=atm_grid_file)
    println(atm_grid_file)
    return atm_grid_file
end

function download_esa_cci(dest_file)
    ## the data access is a bit cumbersome in this case as the website asks for a name and affiliation which requires some tricks
    # get csrfmiddlewaretoken
    csrf_url = "http://products.esa-icesheets-cci.org/products/details/cci_sec_2021.zip/"
    session  = HTTP.get(csrf_url; cookies=true)
    body     = Gumbo.parsehtml(String(session.body))
    inpts    = eachmatch(Selector("input"), body.root)
    function get_token(inpts)
        for inp in inpts
            if !haskey(inp.attributes, "name")
                continue
            end
            if getattr(inp, "name") == "csrfmiddlewaretoken"
                return getattr(inp, "value")
            end
        end
        @error "No token found."
    end
    token = get_token(inpts)

    # send post
    username    = "Peeves"
    affiliation = "Hogwarts"
    credentials = Dict(
        "username" => username,
        "affiliation" => affiliation,
        "csrfmiddlewaretoken" => token
    )
    login_url = "http://products.esa-icesheets-cci.org/register_or_login_no_redirect/"
    HTTP.post(login_url, body=credentials, cookies=true)

    # download
    file_url = "http://products.esa-icesheets-cci.org/products/download/cci_sec_2021.zip"
    out = HTTP.download(file_url, cookies=true)

    # read file out of zip folder
    r            = ZipFile.Reader(out)
    fnames       = [r.files[i].name for i in eachindex(r.files)]
    fi           = findfirst(startswith.(fnames, "Release/CCI_GrIS_RA") .&& endswith.(fnames, "nc"))
    file_to_read = r.files[fi]
    write(dest_file, read(file_to_read, String));
    close(r)
    rm(out)
    return
end

function create_dhdt_grid(gr::Int; startyr::Int, endyr::Int)
    dhdt_dir       = "data/dhdt/"
    download_file  = joinpath(dhdt_dir, "CCI_GrIS_RA_SEC_5km_Vers3.0_2021-08-09.nc")
    get_filename(starty, endy) = splitext(download_file)[1]*"_g$(gr)_$(starty)-$(endy).nc"
    if isfile(get_filename(startyr, endyr))
        return get_filename(startyr, endyr), endyr-startyr
    end
    # download
    if !isfile(download_file)
        println("Downloading dhdt data...")
        download_esa_cci(download_file)
        @assert isfile(download_file)
    end

    # extract cumulative elevation change over certain years
    println("Calculating cumulative elevation change...")
    m0           = NCDataset(download_file)["SEC"]
    t            = NCDataset(download_file)["time"] # center years, elevation change values are moving averages over 4yr time windows
    ti1          = findmin(abs.(Year.(t) .- Year(startyr)))[2]
    actual_start = Year(t[ti1]).value
    if startyr != actual_start @warn "Given start year $startyr was not available, starting at the nearest year $actual_start instead." end
    tin          = findmin(abs.(Year.(t) .- Year(endyr)))[2]
    actual_end   = Year(t[tin]).value
    if endyr != actual_end @warn "Given end year $endyr was not available, ending at the nearest year $actual_end instead." end

    # sum up annual elevation change and save
    n_years = actual_end - actual_start
    msum    = sum(m0[ti1:tin,:,:], dims=1)[1,:,:]
    msum[ismissing.(msum)] .= no_data_value
    tempname = "temp.nc"
    svd_IceSheetDEM.save_netcdf(tempname, download_file, [msum], ["msum"], Dict("msum"=> Dict{String, Any}()))

    # gdalwarp to right grid
    dest = get_filename(actual_start, actual_end)
    gdalwarp(tempname; srcnodata="$no_data_value", gr, dest)[:,end:-1:1]
    rm(tempname)
    return dest, n_years
end

# python functions
function __init__()
    py"""
    import xdem
    import geoutils as gu
    import pandas as pd
    import geopandas as gpd
    import numpy as np
    import verde as vd
    import matplotlib.pyplot as plt

    def point_interp(fname_ref, fname_ref_geoid, fname_atm, fname_out):
        ref_DEM       = xdem.DEM(fname_ref)
        ref_DEM_geoid = xdem.DEM(fname_ref_geoid)
        # extract ATM points
        df       = pd.read_csv(fname_atm)
        geometry = gpd.points_from_xy(df.lon, df.lat, crs="WGS84")
        g        = geometry.to_crs(ref_DEM.crs)
        # interpolate
        ref_pts       = ref_DEM.interp_points(pts=list(zip(g.x, g.y)), prefilter=False, order=2)
        ref_pts_geoid = ref_DEM_geoid.interp_points(pts=list(zip(g.x, g.y)), prefilter=False, order=2)
        # save
        # note: "h_ref" below is referenced to geoid !!
        # (easier to remove geoid from grimp to do the grimp-atm difference rather than add it to atm, because already on the same grid;
        # however, for destandardization we need the geoid-referenced elevation of grimp)
        ds_save = pd.DataFrame({"x": g.x, "y": g.y, "h_ref": ref_pts_geoid, "dh": (ref_pts-df.z)})
        ds_save.to_csv(fname_out, index=False)

    def block_reduce(fname_in, fname_out, spacing):
        df                      = pd.read_csv(fname_in)
        reducer                 = vd.BlockReduce(reduction=np.median, spacing=spacing)
        coordinates, dh_reduced = reducer.filter((df.x, df.y), df.dh)
        coordinates, h_reduced  = reducer.filter((df.x, df.y), df.h_ref)
        xn, yn                  = coordinates
        df_reduced              = pd.DataFrame({"x":xn, "y":yn, "h_ref":h_reduced, "dh":dh_reduced})
        df_reduced.to_csv(fname_out, index=False)

    def coregistration(reference_file, dem_file_not_aligned, dest_file_aligned, outline_shp_file, fig_name):
        dem_not_aligned = xdem.DEM(dem_file_not_aligned)
        reference_dem   = xdem.DEM(reference_file)
        # calculate dh before co-registration
        diff_before     = reference_dem - dem_not_aligned
        # Create a mask of stable terrain, removing outliers outside 3 NMAD
        glacier_outlines = gu.Vector(outline_shp_file)
        mask_noglacier   = ~glacier_outlines.create_mask(reference_dem)
        mask_nooutliers  = np.abs(diff_before - np.nanmedian(diff_before)) < 3 * xdem.spatialstats.nmad(diff_before)
        # Create inlier mask
        inlier_mask      = mask_noglacier & mask_nooutliers
        # co-register
        nuth_kaab = xdem.coreg.NuthKaab()
        print("Doing Nuth and Kaab co-registration...")
        nuth_kaab.fit(reference_dem, dem_not_aligned, inlier_mask)
        print(nuth_kaab._meta)
        aligned_dem = nuth_kaab.apply(dem_not_aligned)
        # calculate dh after co-registration
        diff_after = reference_dem - aligned_dem
        # make a zoomed-in plot to show the difference
        plt.figure(figsize=(14,7))
        plt.subplot(1,2,1)
        plt.imshow(diff_before.data[1600:2150,3500:4300], cmap="coolwarm_r", vmin=-50, vmax=50, interpolation="None")
        plt.colorbar()
        plt.title("before co-registration")
        plt.subplot(1,2,2)
        plt.imshow(diff_after.data[1600:2150,3500:4300], cmap="coolwarm_r", vmin=-50, vmax=50, interpolation="None")
        plt.colorbar()
        plt.title("after co-registration")
        plt.savefig(fig_name)
        # save aligned dem
        dem_xa = aligned_dem.to_xarray("surface")
        dem_xa.to_netcdf(dest_file_aligned)
    """
end
py_point_interp(fname_ref, fname_ref_geoid, fname_atm, fname_out) = py"point_interp"(fname_ref, fname_ref_geoid, fname_atm, fname_out)
py_block_reduce(fname_in, fname_out, spacing) = py"block_reduce"(fname_in, fname_out, spacing)
py_coregistration(reference_file, dem_file_not_aligned, dest_file_aligned, outline_shp_file, fig_name) = py"coregistration"(reference_file, dem_file_not_aligned, dest_file_aligned, outline_shp_file, fig_name)
