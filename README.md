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
- `DEM_projection_GDAL.sh`: bash script to reproject a DEM on the model grid

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
  --λ, --lambda Λ      regularization parameter for the least squares
                       problem (type: Float64, default: 0.0)
  --r, --trunc R       number of modes that should be kept in the SVD
                       truncation (type: Int64)
  --res RES            resolution in m, currently available at '1200'
                       or '1800' (default: "1200")
  --filepath FILEPATH  folder where the training data 'usurf_*' is
                       stored (default: "data/")
  --save               save the output in an .nc file (option takes no
                       argument)
  -h, --help           show this help message and exit
```
So without specifying any arguments it will be run without regularization (λ=0), doing the full SVD (r=nothing) and a resolution of 1200 m (res="1200"). It will also assume that all data is saved in a folder `data/` and it won't save any output. To change that, run it for instance with:
```
$ julia --project SVD_lsqfit.jl --r "100" --res "1800" --save
```


## Run the script interactively in the REPL:
```
$ julia --project
julia> include("SVD_lsqfit.jl")
```
To change the command line arguments from within the REPL, modify the ARGS variable within Julia:
```
$ julia --project
julia> ARGS = ["--r", "100", "--save", "--λ", "34"] # don't forget the "--"!
julia> include("SVD_lsqfit.jl")

# check which arguments it ran with
julia> parse_commandline(ARGS)
Dict{String, Any} with 5 entries:
  "λ"        => 34.0
  "res"      => "1200"
  "r"        => 100
  "save"     => true
  "filepath" => "data/"
  ```
