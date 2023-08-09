# svd_IceSheetDEM

[![CI](https://github.com/pohlan/svd_IceSheetDEM/actions/workflows/CI.yml/badge.svg)](https://github.com/pohlan/svd_IceSheetDEM/actions/workflows/CI.yml)

## Installation / setup
In the shell:
```
$ git clone git@github.com:pohlan/svd_IceSheetDEM.git
$ cd svd_IceSheetDEM
```
Next, download all the required packages either directly in the shell over
```
$ julia --project -e 'import Pkg; Pkg.instantiate()'
```
or in the Julia command line
```
$ julia --project      # --project directly activates the environment
julia> ]               # takes you to the command line of the package manager
(svd_IceSheetDEM) pkg> instantiate
```
## Requirements

### 1) Local files
All training data and a shape file outlining the ice must be available locally.
### 2) Earthdata login
Downloading the data requires earthdata credentials. To provide this information, one can create a file called .netcr in the home folder containing the following:
```
machine urs.earthdata.nasa.gov
login myusername
password abcd123
```
replacing `myusername` with the actual username and `abcd123` with the password.


## Files
- `main.jl`: Downloads the required data and performs preprocessing steps (only the first time). Then it solves the least square fit problem and saves a netcdf in the same format as a bedmachine file.
- `SVD_lsqfit.tex`: tex file deriving the least squares solution
- `params_vs_error.jl`: Plots the mean absolute error of the reconstruction with respect to the input data against different values of a certain parameter (L-curve approach). The parameter can either be the regularization parameter, the number of ensembles in the training data or the truncation of the SVD. Note that this is not the best approach since it is compared to the data that already went into the least square solution, but unfortunately there is no other validation data available.
- `data_download_processing/`: retrieving ATM and bamber data, maybe needed at a later point

## Run the script from the shell
The `main.jl` can be run directly from the shell with the following command line arguments:

```
$ juliap main.jl --help

usage: main.jl [--λ Λ] [--r R] --training_data [TRAINING_DATA...]
               [--imbie_shp_file IMBIE_SHP_FILE] [-h]

optional arguments:
  --λ, --lambda Λ       regularization parameter for the least squares
                        problem (type: Float64, default: 100000.0)
  --r R                 truncation of SVD, default is a full SVD
                        (type: Int64, default: 5076944270305263616)
  --training_data [TRAINING_DATA...]
                        training files, e.g. train_folder/usurf*.nc
  --imbie_shp_file IMBIE_SHP_FILE
                        shape file outlining the ice
  -h, --help            show this help message and exit
```
By default it will thus run with regularization λ=1e5 and a full SVD (which happens when r is larger than the matrix size). The training data and imbie shape files always have to be provided locally and their paths have to be passed as a commandline argument, e.g.:
```
$ julia --project main.jl --training_data training_data/usurf_*.nc --imbie_shp_file imbie_data/gris-outline-imbie-1980_updated.shp
```


## Run the script interactively in the REPL

To change the command line arguments from within the REPL, uncomment the following section in `main.jl` and modify it for the desired parameters:

https://github.com/pohlan/svd_IceSheetDEM/blob/4c36cb9a2046b1a6dfc8f6a7891a7da86403b645/main.jl#L7-L12

Then run

```
$ julia --project
julia> include("main.jl")
```
