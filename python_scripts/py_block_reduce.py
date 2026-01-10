import verde
import pandas as pd
import argparse
import numpy as np

parser = argparse.ArgumentParser()
parser.add_argument("--fname_in",  type=str, help="Input point data")
parser.add_argument("--fname_out", type=str, help="Output point data, block-averaged")
parser.add_argument("--spacing",   type=float, help="LLength over which to average")

# retrieve command line inputs
args = parser.parse_args()
fname_in  = args.fname_in
fname_out = args.fname_out
spacing   = args.spacing

# block reduction
df                      = pd.read_csv(fname_in)
reducer                 = verde.BlockReduce(reduction=np.median, spacing=spacing)
coordinates, dh_reduced = reducer.filter((df.x, df.y), df.dh)
coordinates, h_reduced  = reducer.filter((df.x, df.y), df.h_ref)
xn, yn                  = coordinates
df_reduced              = pd.DataFrame()
df_reduced["x"]         = xn
df_reduced["y"]         = yn
df_reduced["h_ref"]     = h_reduced
df_reduced["dh"]        = dh_reduced
df_reduced.to_csv(fname_out, index=False)
