import xdem
import geoutils as gu
import rioxarray as rio
import pandas as pd
import geopandas as gpd
import argparse

parser = argparse.ArgumentParser()
parser.add_argument("--reference_file",    type=str, help="Raster to coregister to")
parser.add_argument("--point_data_file",   type=str, help="Point data to be coregistered")
parser.add_argument("--outline_shp_file",  type=str, help="Shape file of ice sheet outline")
parser.add_argument("--dest_file_aligned", type=str, help="Output file of coregistered point data")

# retrieve command line inputs
args = parser.parse_args()
reference_file    = args.reference_file
point_data_file   = args.point_data_file
outline_shp_file  = args.outline_shp_file
dest_file_aligned = args.dest_file_aligned

# read reference DEM
ds            = rio.open_rasterio(reference_file)
reference_dem = xdem.DEM.from_xarray(ds)
# read point data
df            = pd.read_csv(point_data_file)
geometry      = gpd.points_from_xy(df.lon, df.lat, crs="WGS84")
g             = geometry.to_crs(reference_dem.crs)
# GeoDataFrame
gdf      = gpd.GeoDataFrame(geometry=g, crs=reference_dem.crs)
gdf["z"] = df.z
gdf["E"] = g.x   # xdem is expecting columns named "E" and "N"
gdf["N"] = g.y
# stable terrain mask (no outlier removal here, in contrast to raster version)
glacier_outlines = gu.Vector(outline_shp_file)
mask_noglacier   = ~glacier_outlines.create_mask(reference_dem)
# coregistration
print("Doing Nuth and Kaab fit...")
nuth_kaab        = xdem.coreg.NuthKaab()
nuth_kaab.fit(reference_dem, gdf, mask_noglacier)
gdf_aligned      = nuth_kaab.apply(gdf)
print(nuth_kaab.meta)
# save
df_save = pd.DataFrame()
df_save["z"] = gdf_aligned.z
df_save["x"] = gdf_aligned.geometry.x
df_save["y"] = gdf_aligned.geometry.y
df_save.to_csv(dest_file_aligned, index=False)
