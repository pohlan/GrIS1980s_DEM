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
    get_options(; grd, cut_shp="", x_min=-678650, y_min=-3371600, x_max=905350, y_max=- 635600)
Returns a vector of Strings that can be used as an input for gdalwarp
"""
function gdalwarp_options(;grd, cut_shp = "", srcnodata="", resampling="bilinear")
    lims = DomainLimits()
    options = ["-overwrite",
               "-t_srs", lims.crs,
               "-r", resampling,
               "-co", "FORMAT=NC4",
               "-co", "COMPRESS=DEFLATE",
               "-co", "ZLEVEL=2",
               "-dstnodata", string(no_data_value),
               "-te", string(lims.x_min), string(lims.y_min), string(lims.x_max), string(lims.y_max),
               "-tr", "$grd", "$grd",
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
    gdalwarp(path; grd, kwargs...)

taken from https://discourse.julialang.org/t/help-request-using-archgdal-gdalwarp-for-resampling-a-raster-dataset/47039

# Example
```
out = gdalwarp("data/input.nc"; grd=1200, dest="data/output.nc")
```
"""
function gdalwarp(path::String; grd::Real, cut_shp="", srcnodata="", resampling="bilinear", kwargs...)
    ds = AG.read(path) do source
        AG.gdalwarp([source], gdalwarp_options(;grd, cut_shp, srcnodata, resampling); kwargs...) do warped
           band = AG.getband(warped, 1)
           AG.read(band)
       end
    end
    return ds
end
function gdalwarp(path::Vector{String}; grd::Real, cut_shp="", srcnodata="", resampling="bilinear", kwargs...)
    source = AG.read.(path)
    ds = AG.gdalwarp(source, gdalwarp_options(;grd, srcnodata, resampling); kwargs...) do warped1
        if !isempty(cut_shp)
            AG.gdalwarp([warped1], gdalwarp_options(;grd, cut_shp, srcnodata, resampling); kwargs...) do warped2
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

function gdalgrid(filenm::String; grd::Real, kwargs...)
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
    rad     = 2*grd
    options = ["-a", "invdist:radius=$rad:min_points=1:nodata=$(string(no_data_value))",
               "-co", "FORMAT=NC4",
               "-co", "COMPRESS=DEFLATE",
               "-co", "ZLEVEL=2",
               "-txe", string(lims.x_min), string(lims.x_max),
               "-tye", string(lims.y_min), string(lims.y_max),
               "-tr", "$grd", "$grd"]
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
    grd = template["x"][2] - template["x"][1]
    ds.attrib["spacing in m"] = "$grd"

    close(ds)
end

function create_bedmachine_grid(grd)
    bedmachine_path = "data/bedmachine/"
    dest_file = joinpath(bedmachine_path, "bedmachine_g$(grd).nc")
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
        new = gdalwarp("NETCDF:"*bedmachine_original*":"*fn;  grd, srcnodata=string(attr_template.attrib["no_data"]), resampling, dest=spatial_template_file)[:,end:-1:1]
        push!(fieldvals, new)  # flip y-axis due to behaviour of ArchGDAL
    end

    # save as netcdf
    attributes = get_attr(attr_template, layernames)
    save_netcdf(dest_file, spatial_template_file, fieldvals, layernames, attributes)
    rm(spatial_template_file)
    return bedmachine_original, dest_file
end

function create_grimpv2(grd_target, grd_coreg, bedmachine_original; kw="")
    data_path = joinpath("data","grimpv2")
    get_dest_file(grd_target) = joinpath(data_path, "grimpv2_geoid_corrected_g$(Int(grd_target)).nc")
    dest_grd_coreg_file       = get_dest_file(grd_coreg)
    dest_grd_target_file      = get_dest_file(grd_target)
    merged_grd_coreg          = joinpath(data_path,"merged_grimpv2_g$(Int(grd_coreg)).nc")

    if isfile(dest_grd_target_file) && isfile(dest_grd_coreg_file)
        return merged_grd_coreg, dest_grd_coreg_file, dest_grd_target_file
    elseif !isfile(dest_grd_coreg_file)

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
        println("Using gdalwarp to merge the grimpv2 mosaics into one DEM, taking a while..")
        grimp_g150  = gdalwarp(grimp_files; grd=grd_coreg, srcnodata="0.0", dest=merged_grd_coreg)
        grimp_g150  = grimp_g150[:,end:-1:1]

        # apply geoid correction
        mask  = gdalwarp("NETCDF:"*bedmachine_original*":mask";  grd=grd_coreg)[:,end:-1:1]
        geoid = gdalwarp("NETCDF:"*bedmachine_original*":geoid";  grd=grd_coreg)[:,end:-1:1]
        idx   = findall(grimp_g150 .!= no_data_value .&& .!ismissing.(geoid) .&& mask .!= 0 .&& mask .!= 4 .&& .!ismissing.(mask))
        grimp_highres = zeros(size(grimp_g150))
        grimp_highres[idx] .= grimp_g150[idx]
        # overwrite merged DEM, here only relevant cells have a value
        sample_path = merged_grd_coreg
        grimp_highres[grimp_highres .<= 0.0] .= no_data_value
        save_netcdf(merged_grd_coreg, sample_path, [grimp_highres], ["surface"], Dict("surface"=>Dict{String,Any}()))
        # correct for geoid and save in different file
        grimp_highres[idx] .-= geoid[idx]
        grimp_highres[grimp_highres .<= 0.0] .= no_data_value
        save_netcdf(dest_grd_coreg_file, sample_path, [grimp_highres], ["surface"], Dict("surface"=>Dict{String,Any}()))
    end

    # gdalwarp to desired grid
    gdalwarp(dest_grd_coreg_file; grd=grd_target, srcnodata=string(no_data_value), dest=dest_grd_target_file)

    return merged_grd_coreg, dest_grd_coreg_file, dest_grd_target_file
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

function create_aerodem(grd_target, grd_coreg, outline_shp_file, bedmachine_original, ref_coreg_file_geoid; kw="")
    aerodem_path              = joinpath("data","aerodem")
    get_aero_file(grd_target) = joinpath(aerodem_path, "aerodem_rm-filtered_geoid-corr_g$(Int(grd_target)).nc")
    # get_rm_file(grd_target)   = aerodem_path * "rm_g$(grd_target).nc"
    aerodem_highres_file      = get_aero_file(grd_coreg)
    aerodem_gr_file           = get_aero_file(grd_target)
    # rm_g150_file      = get_rm_file(grd_coreg)
    # rm_gr_file        = get_rm_file(grd_target)

    if isfile(aerodem_gr_file) && isfile(aerodem_highres_file) #&& isfile(rm_gr_file)
        return aerodem_highres_file, aerodem_gr_file #, rm_gr_file
    elseif !isfile(aerodem_highres_file)

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
        aero_rm_filt     = gdalwarp(aerodem_files; grd=grd_coreg, srcnodata=string(0.0), dest=merged_aero_dest)
        rel_mask         = gdalwarp(     rm_files; grd=grd_coreg, srcnodata=string(0.0), dest=merged_rm_dest)

        # filter for observations where reliability value is low
        aero_rm_filt[rel_mask .< 40] .= no_data_value  # only keep values with reliability of at least xx

        # apply geoid correction
        geoid = gdalwarp("NETCDF:"*bedmachine_original*":geoid";  grd=grd_coreg)
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
        not_aligned_file = splitext(aerodem_highres_file)[1]*"_not_aligned"*splitext(aerodem_highres_file)[2]
        save_netcdf(not_aligned_file, sample_path, [aero_rm_filt[:,end:-1:1]], [layername], attributes)
        # save_netcdf(rm_g150_file, sample_path, [Float32.(rel_mask[:,end:-1:1])], ["reliability mask"], Dict("reliability mask" => Dict{String, Any}()))

        # co-registration
        fig_name = joinpath(fig_path, "coregistration_before_after.png")
        py_coreg_raster(ref_coreg_file_geoid, not_aligned_file, aerodem_highres_file, outline_shp_file, fig_name) # note: output is upside down but will be reversed in gdalwarp below

        rm(merged_aero_dest)
        rm(merged_rm_dest)
        rm(not_aligned_file)
    end

    # gdalwarp to desired grid
    gdalwarp(aerodem_highres_file; grd=grd_target, cut_shp=outline_shp_file, srcnodata=string(no_data_value), dest=aerodem_gr_file)
    # gdalwarp(rm_g150_file; grd=grd_target, srcnodata=string(no_data_value), dest=rm_gr_file)

    return aerodem_highres_file, aerodem_gr_file
end

function create_outline_mask(grd, outline_shp_file, sample_file)
    outline_path = "data/gris-imbie-1980/"
    outline_mask_file = joinpath(outline_path, "outline_mask_g$(grd).nc")
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
    gdalwarp(fname_ones; grd, cut_shp=outline_shp_file, dest=outline_mask_file)
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

function get_atm_dh_file(ref_coreg_file_ellips, ref_coreg_file_geoid, outline_shp_file, blockspacing; kw="")
    # get raw data
    raw_atm_data_file = get_atm_raw_file(kw)

    atm_dh_dest_file = joinpath(dirname(raw_atm_data_file), "dh_ref_minus_atm_blocksize_$(Int(blockspacing))m.csv")
    # if file exists already, do nothing
    if isfile(atm_dh_dest_file)
        return atm_dh_dest_file
    end

    # need a file with only one band for interp_points operation in xdem python;
    # plus need to reference to ellipsoid to be consistent with atm before differencing the two

    # co-registration
    atm_data_aligned  = splitext(raw_atm_data_file)[1]*"_aligned"*splitext(raw_atm_data_file)[2]
    py_coreg_points(ref_coreg_file_ellips, raw_atm_data_file, atm_data_aligned, outline_shp_file)

    # interpolate GrIMP DEM on ATM points and calculate difference
    dh_interpolated_file = joinpath(dirname(raw_atm_data_file), "atm_dh_interpolated.csv")
    py_point_interp(ref_coreg_file_ellips, ref_coreg_file_geoid, atm_data_aligned, dh_interpolated_file)

    # block reduce with python package verde; average data that is heavily oversampled in direction of flight
    println("Doing blockreduce...")
    py_block_reduce(dh_interpolated_file, atm_dh_dest_file, blockspacing)

    # delete all temporary files
    rm(atm_data_aligned)
    rm(dh_interpolated_file)

    return atm_dh_dest_file
end

function create_reconstructed_bedmachine(rec_file, dest)
    # load reconstruction and determine grid size
    surfaceDEM = ncread(rec_file, "surface")
    x          = ncread(rec_file, "x")
    grd         = x[2] - x[1]

    # load bedmachine
    _, bedmachine_file = create_bedmachine_grid(grd)
    bedDEM             = ncread(bedmachine_file, "bed")
    bedm_mask          = ncread(bedmachine_file, "mask")
    ice_mask           = (surfaceDEM .!= no_data_value) .&& (surfaceDEM .> bedDEM)

    # calculate floating mask
    ρw            = 1030
    ρi            = 917
    Pw            = - ρw * bedDEM
    Pi            = ρi * (surfaceDEM - bedDEM)
    floating_mask = (Pw .> Pi) .&& (ice_mask)  # floating where water pressure > ice pressure at the bed

    # calculate mask
    new_mask = ones(eltype(bedm_mask), size(bedm_mask))
    new_mask[bedm_mask .== 0] .= 0 # ocean
    new_mask[ice_mask] .= 2
    new_mask[floating_mask]   .= 3

    # make sure DEM is zero everywhere on the ocean and equal to bed DEM on bedrock
    surfaceDEM[new_mask .== 0] .= no_data_value
    surfaceDEM[new_mask .== 1] .= bedDEM[new_mask .== 1]

    # calculate ice thickness
    h_ice = zeros(size(surfaceDEM)) .+ no_data_value
    h_ice[ice_mask] .= surfaceDEM[ice_mask] - bedDEM[ice_mask]
    h_ice[floating_mask] .= surfaceDEM[floating_mask] ./  (1-ρi/ρw)

    # save to netcdf file
    layers      = [surfaceDEM, bedDEM, h_ice, new_mask]
    layernames  = ["surface", "bed", "thickness", "mask"]
    template    = NCDataset(bedmachine_file)
    attributes  = get_attr(template, layernames)
    # overwrite some attributes
    sources_rec = Dict("surface"   => "svd reconstruction",
                       "bed"       => "Bedmachine-v5: Morlighem et al. (2022). IceBridge BedMachine Greenland, Version 5. Boulder, Colorado USA. NASA National Snow and Ice Data Center Distributed Active Archive Center. https://doi.org/10.5067/GMEVBWFLWA7X; projected on new grid with gdalwarp",
                       "thickness" => "computed from surface and bed",
                       "mask"      => "bedrock from Morlighem et al. (2022); ice, floating and ocean computed from surface and bed elevation"
                       )
    for l in layernames
        attributes[l]["source"] = sources_rec[l]
    end
    attributes["mask"]["long_name"] = "mask (0 = ocean, 1 = ice-free land, 2 = grounded ice, 3 = floating ice)"
    save_netcdf(dest, bedmachine_file, layers, layernames, attributes)

    return dest
end

# python functions
const xdem = Ref{Py}()
const np   = Ref{Py}()
const pd   = Ref{Py}()
const gpd  = Ref{Py}()
function __init__()
    xdem[] = pyimport("xdem")
    np[] = pyimport("numpy")
    pd[] = pyimport("pandas")
    gpd[] = pyimport("geopandas")
end

function py_point_interp(fname_ref, fname_ref_geoid, fname_atm, fname_out)
    ref_DEM       = xdem[].DEM(fname_ref)
    ref_DEM_geoid = xdem[].DEM(fname_ref_geoid)
    # extract ATM points
    df       = pd[].read_csv(fname_atm)
    geometry = gpd[].points_from_xy(df.x, df.y, crs=ref_DEM.crs)
    g        = geometry.to_crs(ref_DEM.crs)
    # interpolate
    ref_pts       = ref_DEM.interp_points(points=(g.x, g.y), prefilter=false)
    ref_pts_geoid = ref_DEM_geoid.interp_points(points=(g.x, g.y), prefilter=false)
    # save
    # note: "h_ref" below is referenced to geoid !!
    # (easier to remove geoid from grimp to do the grimp-atm difference rather than add it to atm, because already on the same grid;
    # however, for destandardization we need the geoid-referenced elevation of grimp)
    ds_save = pd[].DataFrame()
    ds_save["x"]     = g.x
    ds_save["y"]     = g.y
    ds_save["h_ref"] = ref_pts_geoid
    ds_save["dh"]    = ref_pts-df.z
    ds_save.to_csv(fname_out, index=false)
end

function py_block_reduce(fname_in, fname_out, spacing)
    df                      = pd[].read_csv(fname_in)
    reducer                 = pyimport("verde").BlockReduce(reduction=np[].median, spacing=spacing)
    coordinates, dh_reduced = reducer.filter((df.x, df.y), df.dh)
    coordinates, h_reduced  = reducer.filter((df.x, df.y), df.h_ref)
    xn, yn                  = coordinates
    df_reduced              = pd[].DataFrame()
    df_reduced["x"]         = xn
    df_reduced["y"]         = yn
    df_reduced["h_ref"]     = h_reduced
    df_reduced["dh"]        = dh_reduced
    df_reduced.to_csv(fname_out, index=false)
end

function py_coreg_raster(reference_file, dem_file_not_aligned, dest_file_aligned, outline_shp_file, fig_name)
    dem_not_aligned = xdem[].DEM(dem_file_not_aligned)
    reference_dem   = xdem[].DEM(reference_file)
    # calculate dh before co-registration
    diff_before     = reference_dem - dem_not_aligned
    # Create a mask of stable terrain, removing outliers outside 3 NMAD
    glacier_outlines = pyimport("geoutils").Vector(outline_shp_file)
    mask_noglacier   = ~glacier_outlines.create_mask(reference_dem)
    mask_nooutliers  = np[].abs(diff_before - np[].nanmedian(diff_before)) < 3 * xdem[].spatialstats.nmad(diff_before)
    # Create inlier mask
    inlier_mask      = mask_noglacier & mask_nooutliers
    # co-register
    nuth_kaab = xdem[].coreg.NuthKaab()
    println("Doing Nuth and Kaab co-registration...")
    nuth_kaab.fit(reference_dem, dem_not_aligned, inlier_mask)
    println(nuth_kaab._meta)
    aligned_dem = nuth_kaab.apply(dem_not_aligned)
    # calculate dh after co-registration
    diff_after = reference_dem - aligned_dem
    # make a zoomed-in plot to show the difference
    diff_before_jl = PyArray(diff_before.data); diff_before_jl[PyArray(diff_before.data.mask)] .= NaN
    diff_after_jl  = PyArray(diff_after.data);    diff_after_jl[PyArray(diff_after.data.mask)] .= NaN
    p1 = heatmap(diff_before_jl[1600:2150,3500:4300], cmap=:coolwarm, clims=(-50,50), title="before co-registration", aspect_ratio=1, size=(700,900))
    p2 = heatmap(diff_after_jl[1600:2150,3500:4300], cmap=:coolwarm, clims=(-50,50), title="after co-registration", aspect_ratio=1, size=(700,900))
    plot(p1,p2,size=(1500,900))
    savefig(fig_name)
    # save aligned dem
    dem_xa = aligned_dem.to_xarray("surface")
    dem_xa.to_netcdf(dest_file_aligned)
end

# without calculating and plotting diff_before and diff_after
function py_coreg_points(reference_file, point_data_file, dest_file_aligned, outline_shp_file)
    # read reference DEM
    ds            = pyimport("rioxarray").open_rasterio(reference_file)
    reference_dem = xdem[].DEM.from_xarray(ds)
    # read point data
    df            = pd[].read_csv(point_data_file)
    geometry      = gpd[].points_from_xy(df.lon, df.lat, crs="WGS84")
    g             = geometry.to_crs(reference_dem.crs)
    # GeoDataFrame
    gdf      = gpd[].GeoDataFrame(geometry=g, crs=reference_dem.crs)
    gdf["z"] = df.z
    gdf["E"] = g.x   # xdem is expecting columns named "E" and "N"
    gdf["N"] = g.y
    # stable terrain mask (no outlier removal here, in contrast to raster version)
    glacier_outlines = pyimport("geoutils").Vector(outline_shp_file)
    mask_noglacier   = ~glacier_outlines.create_mask(reference_dem)
    # coregistration
    println("Doing Nuth and Kaab fit...")
    nuth_kaab        = xdem[].coreg.NuthKaab()
    nuth_kaab.fit(reference_dem, gdf, mask_noglacier)
    gdf_aligned      = nuth_kaab.apply(gdf)
    println(nuth_kaab.meta)
    # save
    df_save = pd[].DataFrame()
    df_save["z"] = gdf_aligned.z
    df_save["x"] = gdf_aligned.geometry.x
    df_save["y"] = gdf_aligned.geometry.y
    df_save.to_csv(dest_file_aligned, index=false)
    return
end
