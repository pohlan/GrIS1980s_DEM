import xdem
import geoutils as gu
import pandas as pd
import geopandas as gpd
import argparse

parser = argparse.ArgumentParser()
parser.add_argument("--aerodem_file", type=str, help="AeroDEM")
parser.add_argument("--atm_file",     type=str, help="ATM points (ellipsoid)")
parser.add_argument("--output_file",  type=str, help="Output csv file for aerodem_minus_atm")

# retrieve command line inputs
args = parser.parse_args()
aerodem_file = args.aerodem_file
atm_file     = args.atm_file
output_file  = args.output_file

# atm_file   = "data/ATM/dh_ref_minus_atm_blocksize_400m.csv"
# aerodem_file  = "data/aerodem/aerodem_rm-filtered_geoid-corr_g150.nc"
# output_file   = "data/ATM/aero_minus_ATM.csv"

DEM_aero = xdem.DEM(aerodem_file)

# extract ATM points
df       = pd.read_csv(atm_file)
g = gpd.points_from_xy(df.x, df.y, crs=DEM_aero.crs)
# interpolate
aerodem_pts       = DEM_aero.interp_points(points=(g.x, g.y), prefilter=False)
# get h_ATM (not saved directly)
h_atm = df.h_ref - df.dh
# save
ds_save = pd.DataFrame()
ds_save["x"]     = g.x
ds_save["y"]     = g.y
ds_save["dh"]    = aerodem_pts-h_atm
ds_save.to_csv(output_file, index=False)
