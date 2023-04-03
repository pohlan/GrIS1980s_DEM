# svd_IceSheetDEM

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
## Files
- `SVD_lsqfit.jl`: Solves the least square fit problem, see below how to run it
- `SVD_lsqfit.tex`: tex file deriving the least squares solution applied in the .jl script
- `lambda_vs_error.jl`: Plots the mean absolute error of the reconstruction with respect to the input data against different values of the regularization parameter lambda (in order to find a suitable lambda). Not the best approach since it is compared to the data that already went into the least square solution, but unfortunately there is no other validation data available.
- `data_download_processing/`: files to rejproject DEMs to model grid, apply geoid corrections etc. as well as downloading routine for ATM data

## Run the script from the shell
The `SVD_lsqfit.jl` can be run directly from the shell with
```
$ julia --project SVD_lsqfit.jl
```
It comes with the following options for command line arguments:

```
$ julia --project SVD_lsqfit.jl --help

usage: SVD_lsqfit.jl [--λ Λ] [--r R] [--res RES] [--filepath FILEPATH]
                     [--save] [-h]

optional arguments:
  --λ, --lambda Λ       regularization parameter for the least squares
                        problem (type: Float64, default: 100000.0)
  --res RES             resolution in m, currently available at '1200'
                        or '1800' (default: "1200")
  --train_folder TRAIN_FOLDER
                        folder where the training data 'usurf_*' is
                        stored (default: "data/")
  --obs OBS             file of observations that the SVD is fitted to
                        (default:
                        "data/aerodem_g1200m_geoid_corrected_1978_1987_mean.nc")
  --obs_band_name OBS_BAND_NAME
                        name of the surface elevation band in the
                        netcdf file (default: "surface_altitude")
  --save                save the output in an .nc file (option takes
                        no argument)
  -h, --help            show this help message and exit
```
So without specifying any arguments it will be run without regularization (λ=0), doing the full SVD (r=nothing) and a resolution of 1200 m (res="1200"). It will also assume that all data is saved in a folder `data/` and it won't save any output. To change that, run it for instance with:
```
$ julia --project SVD_lsqfit.jl --res "1800" --save
```


## Run the script interactively in the REPL:
```
$ julia --project
julia> include("SVD_lsqfit.jl")
```
To change the command line arguments from within the REPL, modify the ARGS variable within Julia:
```
$ julia --project
julia> ARGS = ["--save", "--λ", "500"] # don't forget the "--"!
julia> include("SVD_lsqfit.jl")

# check which arguments it ran with
julia> parse_commandline(ARGS)
Dict{String, Any} with 6 entries:
  "obs"           => "data/aerodem_g1200m_geoid_corrected_1978_1987_mean.nc"
  "obs_band_name" => "surface_altitude"
  "λ"             => 500.0
  "res"           => "1200"
  "save"          => true
  "train_folder"  => "data/"
  ```
