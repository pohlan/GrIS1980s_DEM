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

## Solving the SVD least-square problem
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

To run the script interactively in the REPL:
```
$ julia --project
julia> include("SVD_lsqfit.jl")
```