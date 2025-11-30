using GrIS1980s_DEM, NCDatasets, Meshes, GeoStats, CSV, DataFrames, JLD2, UnPack, ProgressMeter

# set target directories
main_output_dir = joinpath("output","validation")
mkpath(main_output_dir)
fig_dir         = joinpath(main_output_dir, "figures")
mkpath(fig_dir)

# for running the script interactively
# ARGS = [
#         "--shp_file", "data/gris-imbie-1980/gris-outline-imbie-1980_updated.shp",
#         "--grid_size", 600.0]

parsed_args         = parse_commandline(ARGS)
outline_shp_file    = parsed_args["shp_file"]
grd                 = parsed_args["grid_size"]

################################################
# Preprocessing, standardization and variogram #
################################################

csv_preprocessing, jld2_preprocessing = prepare_obs(grd, outline_shp_file)

# get I_no_ocean, (de-)standardization functions and variogram from pre-processing
df_all = CSV.read(csv_preprocessing, DataFrame)
coords_obs = [GrIS1980s_DEM.F.([x,y]) for (x,y) in zip(df_all.x, df_all.y)]
dict   = load(jld2_preprocessing)
@unpack I_no_ocean, idx_aero, gamma, gamma_error, href_file = dict
varg = GrIS1980s_DEM.get_var(gamma)
varg_error = GrIS1980s_DEM.get_var(gamma_error, nVmax=1)
kernel_signal = varg_to_kernel(varg)
kernel_error  = varg_to_kernel(varg_error)

# make geotable
i_atm = findall(df_all.source .== "atm")
geotable_all = GrIS1980s_DEM.make_geotable(df_all.dh_detrend, df_all.x, df_all.y)
geotable_atm = GrIS1980s_DEM.make_geotable(df_all.dh_detrend[i_atm], df_all.x[i_atm], df_all.y[i_atm])

# add measurement uncertainty for AeroDEM / ATM
GrIS1980s_DEM.add_sigma_obs!(df_all)

##############################################
# Standard cross-validation, leave block out #
##############################################

# create sets of training and test data
ℓ         = 2e5
δl        = 1.5e5
blocks    = partition(geotable_atm, BlockPartition(ℓ))
ids_test  = indices(blocks)
ids_train = []
for i in eachindex(ids_test)
    ids_test[i] = i_atm[ids_test[i]]
    x_min, x_max = extrema(first.(coords_obs[ids_test[i]]))
    y_min, y_max = extrema(last.(coords_obs[ids_test[i]]))
    i_big = findall(x_min-δl .<= df_all.x .<= x_max+δl .&& y_min-δl .<= df_all.y .<= y_max+δl)
    push!(ids_train, setdiff(i_big, ids_test[i])) # indices that are in i_big but not in ids_test[i]
end

println("GP cross-validation...")
function evaluate_fun(i_train, i_test)
    m_pred = GrIS1980s_DEM.do_GP(coords_obs[i_train], GrIS1980s_DEM.F.(df_all.dh_detrend[i_train]), coords_obs[i_test], kernel_signal, kernel_error, df_all.sigma_obs[i_train], var=false)
    return m_pred
end
difs = GrIS1980s_DEM.step_through_folds(ids_train, ids_test, evaluate_fun, df_all.dh_detrend, save_distances=true)

# calculate distance to closest observation
dists = nearest_neighb_distance_from_cv(ids_train, ids_test, df_all.x, df_all.y)

# get indices and save
idx  = vcat(ids_test...)
logℓ = round(log10(ℓ),digits=1)
dest = get_cv_file_GP(grd, logℓ)
cv_dict = (; difs, dists, idx, binfield1=df_all.bfield_1[idx], h_ref=df_all.h_ref[idx], grd, method="GP")
jldsave(dest; cv_dict...)

# # make map of distances to closest observation (for uncertainty estimation)
# ir_sim      = setdiff(I_no_ocean, idx_aero)
# x           = NCDataset(href_file)["x"][:]; nx = length(x)
# y           = NCDataset(href_file)["y"][:]
# x_sim       = x[get_ix.(ir_sim, nx)]
# y_sim       = y[get_iy.(ir_sim, nx)]
# coords_sim  = [GrIS1980s_DEM.F.([x,y]) for (x,y) in zip(x_sim, y_sim)]

# min_dists = zeros(length(x),length(y)) .+ GrIS1980s_DEM.no_data_value
# min_dists[ir_sim] = nearest_neighb_distance_raster(ir_sim, coords_sim, geotable_all)
# dest_file = get_distance_file()
# GrIS1980s_DEM.save_netcdf(dest_file, href_file, [min_dists], ["distance"], Dict("distance" => Dict{String,Any}()))
