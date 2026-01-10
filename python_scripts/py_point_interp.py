import xdem
import pandas as pd
import geopandas as gpd
import argparse

parser = argparse.ArgumentParser()
parser.add_argument("--ref_file",        type=str, help="Raster to interpolate, on ellipsoid")
parser.add_argument("--ref_file_geoid",  type=str, help="Raster to interpolate, on geoid")
parser.add_argument("--point_data_file", type=str, help="Points to interpolate raster on")
parser.add_argument("--dest_file",       type=str, help="Output file of interpolated point data")

# retrieve command line inputs
args = parser.parse_args()
ref_file        = args.ref_file
ref_file_geoid  = args.ref_file_geoid
point_data_file = args.point_data_file
dest_file       = args.dest_file

# load
ref_DEM       = xdem.DEM(ref_file)
ref_DEM_geoid = xdem.DEM(ref_file_geoid)
# extract ATM points
df       = pd.read_csv(point_data_file)
geometry = gpd.points_from_xy(df.x, df.y, crs=ref_DEM.crs)
g        = geometry.to_crs(ref_DEM.crs)
# interpolate
ref_pts       = ref_DEM.interp_points(points=(g.x, g.y), prefilter=False)
ref_pts_geoid = ref_DEM_geoid.interp_points(points=(g.x, g.y), prefilter=False)
# save
# note: "h_ref" below is referenced to geoid !!
# (easier to remove geoid from grimp to do the grimp-atm difference rather than add it to atm, because already on the same grid;
# however, for destandardization we need the geoid-referenced elevation of grimp)
ds_save = pd.DataFrame()
ds_save["x"]     = g.x
ds_save["y"]     = g.y
ds_save["h_ref"] = ref_pts_geoid
ds_save["dh"]    = ref_pts-df.z
ds_save.to_csv(dest_file, index=False)
