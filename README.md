# svd_IceSheetDEM

## Installation / setup
In the shell:
```
$ git clone git@github.com:pohlan/svd_IceSheetDEM.git
$ cd svd_IceSheetDEM
$ julia --project      # --project directly activates the environment
```
In Julia, download all the required packages via
```
julia> import Pkg; Pkg.instantiate()
```

## Solving the SVD least-square problem
Open `SVD_lsqfit.jl` and adjust the file paths to wherever the data is (definition of `model_files` variable).

Then
```
julia> include("SVD_lsqfit.jl")
```
