import xdem
import geoutils as gu
import numpy as np
from pathlib import Path
import argparse

parser = argparse.ArgumentParser()
parser.add_argument("--reference_file",       type=str, help="Reference DEM to coregister to")
parser.add_argument("--dem_file_not_aligned", type=str, help="Input file of DEM to be coregistered")
parser.add_argument("--dest_file_aligned",    type=str, help="Output file of aligned DEM")
parser.add_argument("--outline_shp_file",     type=str, help="Shape file of ice sheet outline")
parser.add_argument("--fig_name",             type=str, help="File name to save figure before/after")

# retrieve command line inputs
args = parser.parse_args()
reference_file       = args.reference_file
dem_file_not_aligned = args.dem_file_not_aligned
dest_file_aligned    = args.dest_file_aligned
outline_shp_file     = args.outline_shp_file
fig_name = args.fig_name

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
diff_before.plot
# p1 = heatmap(diff_before_jl[1600:2150,3500:4300], cmap=:coolwarm, clims=(-50,50), title="before co-registration", aspect_ratio=1, size=(700,900))
# p2 = heatmap(diff_after_jl[1600:2150,3500:4300], cmap=:coolwarm, clims=(-50,50), title="after co-registration", aspect_ratio=1, size=(700,900))
# plot(p1,p2,size=(1500,900))
# savefig(fig_name)
# save aligned dem
dem_xa = aligned_dem.to_xarray("surface")
dem_xa.to_netcdf(dest_file_aligned)
