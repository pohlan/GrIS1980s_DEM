using svd_IceSheetDEM

using Cascadia, Gumbo, HTTP, DataFrames, Downloads
using DataStructures: OrderedDict

dest_path = "data/aerodem_raw/"
mkpath(dest_path)

"""
copied from https://gist.github.com/scls19fr/9ea2fd021d5dd9a97271da317bff6533
"""
function get_table_from_html(input::AbstractString)
    r = HTTP.request("GET", input)
    strr=String(r.body)
    h = Gumbo.parsehtml(strr)
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

# Download
println("Downloading aerodem tif files, this may take a few minutes...")
url_DEMs          = "https://www.nodc.noaa.gov/archive/arc0088/0145405/1.1/data/0-data/G150AERODEM/DEM/"
url_reliablt_mask = "https://www.nodc.noaa.gov/archive/arc0088/0145405/1.1/data/0-data/G150AERODEM/ReliabilityMask/"
aerodem_files     = []
rm_files          = []
for url in [url_DEMs, url_reliablt_mask]
    df = get_table_from_html(url)
    tif_files = df.Name[endswith.(df.Name, ".tif") .&& .!occursin.("carey", df.Name)]
    missing_files = tif_files[.!isfile.(tif_files)]
    Downloads.download.(url .* missing_files, dest_path .* missing_files)
    if startswith(tif_files[1], "aerodem")
        aerodem_files = tif_files
    elseif startswith(tif_files[1], "rm")
        rm_files      = tif_files
    end
end

# gdalwarp
grid = 150
cut_shp = "data/gris-imbie-1980/gris-outline-imbie-1980.shp"
merged_aero_dest = dest_path*"merged_aerodem_g$grid.nc"
merged_rm_dest   = dest_path*"merged_rm_g$grid.nc"
aero     = gdalwarp(dest_path.*aerodem_files; grid, cut_shp, dest=merged_aero_dest)
rel_mask = gdalwarp(dest_path.*     rm_files; grid, cut_shp, dest=merged_rm_dest)

# filter for observations where reliability value is low
aero_rm_filt = copy(aero)
aero_rm_filt[rel_mask .< 40] .= 0.0  # only keep values with reliability of at least xx

# apply geoid correction
geoid                   = shortread("data/bedm_geoid_g$grid.nc")
idx                     = findall(aero_rm_filt .!= 0.0)
aero_rm_geoid_corr      = zeros(size(aero_rm_filt))
aero_rm_geoid_corr[idx] = aero_rm_filt[idx] - geoid[idx]
aero_rm_geoid_corr[aero_rm_geoid_corr .< 0.0] .= 0.0

# save as netcdf
sample_path = merged_aero_dest
dest        = dest_path*"aerodem_rm-filtered_geoid-corr_g$(grid).nc"
save_netcdf(aero_rm_geoid_corr; dest, sample_path)

# move to different resolution
# grid = 1200
# gdalwarp(dest; grid, srcnodata="0.0", dest=dest_path*"aerodem_rm-filtered_geoid-corr_g$(grid).nc")
