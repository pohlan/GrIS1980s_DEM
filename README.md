# GrIS1980s_DEM


## Installation / setup

This code is mostly written in Julia but relies on a few python packages for pre-processing. If pre-processing is done already, only Julia needs to be set up. Otherwise, follow the instructions below to set up the two independent Julia and Python environments (assumes both Python and Julia are installed):

### 1) Julia
Clone the repo:
```
$ git clone git@github.com:pohlan/GrIS1980s_DEM.git
$ cd GrIS1980s_DEM
```
Next, download all the required packages either directly in the shell over
```
$ julia --project -e 'import Pkg; Pkg.instantiate()'
```
or in the Julia command line
```
$ julia --project      # --project directly activates the environment
julia> ]               # takes you to the command line of the package manager
(GrIS1980s_DEM) pkg> instantiate
```
### 2) Python
First, create a conda enviroment using the environment.yml file
```
conda env create -f environment.yml
```
Get the environment path through
```
$ conda activate GrISenv
$ which python
/home/.../python
```
Copy the path and use it as a command line argument 'GrISenv' when running the Julia script (see below).

**Note**: The GrISenv has to **not** be activated when running the Julia script.


## Additional requirements

### 1) Local files
Two things must be available locally:

- a shape file outlining the ice
- the 'training data' (PISM generated realizations of surface elevation in the study)

### 2) Earthdata login
Downloading the data requires earthdata credentials. To provide this information, one can create a file called `.netcr` in the home folder containing the following:
```
machine urs.earthdata.nasa.gov
login myusername
password abcd123
```
replacing `myusername` with the actual username and `abcd123` with the password.


## Running the scripts

### Scripts to run full analysis from scratch (probably don't want to run everything in one go)
The first script that is run will download all the data and do the pre-processing. Scripts can be run independently except `plot_all.jl` obviously as it relies on the output.
1) `cross_validation_GP.jl`: the 'shortest' but may still take several hours / a full day
2) `cross_validation_SVD.jl`: high memory requirement (>64GB RAM) depending on the amount of 'training data'
3) `GP_reconstruction.jl`: may take several days
4) `SVD_reconstruction.jl`: high memory requirement (>64GB RAM) depending on the amount of 'training data'
4) `plot_scripts/plot_all.jl`: produce figures for the study, all from saved outputs

### Run a script from the shell
The scripts above can be run directly from the shell with the following command line arguments:

```
$ julia --project cross_validation_GP.jl --help

usage: cross_validation_GP.jl [--GrISenv GRISENV] [--λ Λ] [--r R]
                        [--training_data [TRAINING_DATA...]]
                        [--shp_file SHP_FILE]
                        [--use_arpack USE_ARPACK] [--maxn MAXN]
                        [--grid_size GRID_SIZE] [-h]

optional arguments:
  --GrISenv GRISENV     output of `which python`; needed for calling
                        python scripts in co-registration
  --λ, --lambda Λ       regularization parameter for the least squares
                        problem (type: Float64, default: 1.0e7)
  --r R                 truncation of SVD, default is a full SVD
                        (type: Int64, default: 5076944270305263616)
  --training_data [TRAINING_DATA...]
                        model-generated realizations of ice sheet
                        elevation as netcdf files, e.g.
                        train_folder/usurf*.nc
  --shp_file SHP_FILE   shape file outlining the ice sheet
  --use_arpack USE_ARPACK
                        if set to true, the Arpack svd is used instead
                        of the standard LinearAlgebra algorithm;
                        Arpack is iterative and matrix free and thus
                        useful when memory becomes limiting, but it
                        can be slower (type: Bool, default: false)
  --grid_size GRID_SIZE
                        cell size of grid, same in x and y direction;
                        not needed for svd reconstruction where
                        training_data is provided (type: Int64,
                        default: 600)
  -h, --help            show this help message and exit
```
The training data and the shape file outlining the ice sheet have to be provided locally and their paths have to be passed as a commandline argument.

For Gaussian process regression (GP), only the `grid_size` and `shp_file` arguments have to be provided. `GrISenv` has to be provided with the first file that is run and that does the pre-processing.
For example:
```
$ julia --project cross_validation_GP.jl --grid_size 600 --shp_file outline.shp --GrISenv /home/.../envs/GrISenv/bin/python
$ julia --project cross_validation_SVD.jl --training_data training_data/usurf_*.nc --shp_file outline.shp --r 1000 --lambda 1e7
```


## Run the script interactively in the REPL

To prescribe those parameters in interactive mode, every of the above mentioned script has a section that can be uncommented and adjusted according to the desired parameters. For instance for `cross_validation_GP.jl`:

```
# for running the script interactively
# ARGS = [
#         "--shp_file", "data/gris-imbie-1980/gris-outline-imbie-1980_updated.shp",
#         "--grid_size", 600.0]
```

Then run

```
$ julia --project
julia> include("cross_validation_GP.jl")
```
Although, it should be noted that interactive mode is only useful for development as the performance may be reduced.
