using GrIS1980s_DEM
using Test, LinearAlgebra, NetCDF, NCDatasets, CSV, DataFrames, Glob, UnPack, JLD2, GeoStats

#####################################################################################
# This script tests the pre-processing as well as the reconstructions               #
# and cross-validations for both the GP and the SVD.                                #
# This is just a test with a lot less data to check quickly if the whole routine    #
# still does the same thing as it used to.                                          #
#                                                                                   #
# Run the script with e.g. (find the path with 'conda activate' and 'which python') #
# $ conda deactivate                                                                #
# $ julia --project test/runtests.jl --GrISenv /home/.../GrISenv/bin/python         #
#####################################################################################

cd(@__DIR__)  # set working directory to where file is located

const grd = 4000
blockspacing=Int(0.5*grd)
outline_shp_file            = joinpath("testdata", "testshape.shp")

# remove all files with this grid that exist, and any other output in the test folder
for file in glob("**/*g$(grd).*", joinpath(pkgdir(GrIS1980s_DEM), "data"))
    rm(file, recursive=true, force=true)
end
for file in glob("**/*.*", joinpath(pkgdir(GrIS1980s_DEM), "test", "output"))
    rm(file, recursive=true, force=true)
end
rm(joinpath(pkgdir(GrIS1980s_DEM), "data", "ATM", "atm_dh_interpolated.csv"), recursive=true, force=true)
rm(joinpath(pkgdir(GrIS1980s_DEM), "data", "ATM", "ATM_nadir2seg_all_aligned.csv"), recursive=true, force=true)
rm(joinpath(pkgdir(GrIS1980s_DEM), "data", "ATM", "dh_ref_minus_atm_blocksize_$(blockspacing)m.csv"), recursive=true, force=true)

##################
# Pre-processing #
##################

# bedmachine
bedmachine_original, bedmachine_file_gr = GrIS1980s_DEM.create_bedmachine_grid(grd)
# GrIMP v2
ref_coreg_file_ellips, ref_coreg_file_geoid, ref_file   = GrIS1980s_DEM.create_grimpv2(grd, 150, bedmachine_original)
# aerodem
aero_150_file, aero_gr_file = GrIS1980s_DEM.create_aerodem(grd, 150, outline_shp_file, bedmachine_original, ref_coreg_file_geoid)    # only download the DEMs from 1981 while testing to save some space
# imbie mask
outline_mask_file = GrIS1980s_DEM.create_outline_mask(grd, outline_shp_file, aero_150_file)
# atm
atm_dh_file = GrIS1980s_DEM.get_atm_dh_file(ref_coreg_file_ellips, ref_coreg_file_geoid, outline_shp_file, blockspacing)

missmax(x) = maximum(x[.!ismissing.(x)])
missmin(x) = minimum(x[.!ismissing.(x)])
missum(x)  = sum(x[.!ismissing.(x)])

@testset "bedmachine" begin
    # read in
    ds = NCDataset(bedmachine_file_gr)
    @test all(["mask","geoid","bed","surface","thickness"] .∈ (keys(ds),))
    mask  = ds["mask"][:,:]
    geoid = ds["geoid"][:,:]
    bed   = ds["bed"][:,:]
    surface = ds["surface"][:,:]
    # test that the values are correct
    @test missum(mask.==1) .== 20053 && missum(mask.==2) .== 120513 && missum(mask.==3) .== 405 && missum(mask.==4) .== 7570
    @test missmax(geoid) == 64 && missmin(geoid) == 6 && sum(ismissing.(geoid)) == 8208
    @test missmax(bed) ≈ 3010.345f0 && missmin(bed) ≈ -5505.15f0
    @test missmax(surface) ≈ 3233.3943f0 && missmin(surface) ≈ 0.0f0
    # check that none of the layers are "upside down" as happens sometimes with gdalwarp
    ny_half = fld(size(mask,2),2)
    @test missum(mask[:,1:ny_half] .== 4) == 0
    @test missum(geoid[:,1:ny_half]) > missum(geoid[:,(ny_half+1):end])           # more higher values for the geoid in the southern half
    @test missum(bed[:,1:ny_half] .< 0) > missum(bed[:,(ny_half+1):end] .< 0)     # more and deeper ocean bed in the south
    @test missum(surface[:,1:ny_half] .== 0) > missum(surface[:,(ny_half+1):end] .== 0)     # more ocean than in the south
    close(ds)
end

@testset "reference elevation grimp v2" begin
    ds = NCDataset(ref_file)["Band1"][:,:]
    @test missmax(ds) ≈ 3234.319f0
    @test sum(.!ismissing.(ds)) == 133843
    # check that it's not "upside down"
    ny_half = fld(size(ds,2),2)
    @test sum(.!ismissing.(ds[:,(ny_half+1):end])) > sum(.!ismissing.(ds[:,1:ny_half]))   # there should be more data in the northern half
end

@testset "aerodem" begin
    ds150                 = NCDataset(aero_150_file)["surface"][:,:]
    ds150                 = nomissing(ds150,0.0)
    dsgr                  = NCDataset(aero_gr_file )["Band1"][:,:]
    @test missmax(ds150)  ≈ 3639.633f0
    @test sum(.!ismissing.(dsgr)) == 1559
    @test missmax(dsgr) ≈ 2205.294f0
    # check that it's not "upside down"
    ny_half = fld(size(dsgr,2),2)
    sum(.!ismissing.(dsgr[:,(ny_half+1):end])) > sum(.!ismissing.(dsgr[:,1:ny_half]))    # there should be more data in the northern half
end

@testset "outline mask" begin
    outline_mask = NCDataset(outline_mask_file)["Band1"][:,:]
    @test missum(outline_mask .== 1) == sum(.!ismissing.(outline_mask)) == 14627
    # check that it's not "upside down"
    ny_half = fld(size(outline_mask,2),2)
    sum(.!ismissing.(outline_mask[:,(ny_half+1):end])) > sum(.!ismissing.(outline_mask[:,1:ny_half]))  # there should be more data in the northern half
end

@testset "atm" begin
    df_atm = CSV.read(atm_dh_file, DataFrame)
    @test size(df_atm) == (27047,4)
    @test all(extrema(df_atm.x) .≈ (-633871f0, 670411f0))
    @test all(extrema(df_atm.y) .≈ (-3.294190f6, -773827f0))
    @test sum(.!ismissing.(df_atm.dh)) == 26494  # the ones outside of the outline are missing because there is no h_ref
    @test missmax(df_atm.h_ref) ≈ 3235f0
    @test missmin(df_atm.dh) ≈ -1755.70f0 && missmax(df_atm.dh) ≈ 626.52f0
end


####################
# GP INTERPOLATION #
####################

# set the seed, necessary because it randomly selects a portion of aerodem data for variogram
using Random
Random.seed!(1234)

csv_preprocessing, jld2_preprocessing = GrIS1980s_DEM.prepare_obs(grd, outline_shp_file, blockspacing=blockspacing, nbins1=5, nbins2=5, min_n_sample=80)
rec_file = geostats_interpolation(grd, outline_shp_file, csv_preprocessing, jld2_preprocessing, ℓ_block=1e5, δl=1.4e5)

@testset "GP reconstruction" begin
    # surface
    h_GP = NCDataset(rec_file)["surface"][:,:]
    ny_half = fld(size(h_GP,2),2)
    # test that it's not "upside down"
    @test sum(.!ismissing.(h_GP[:,(ny_half+1):end])) == 0
    @test sum(.!ismissing.(h_GP[:,1:ny_half])) == 14623
    # test values
    @test missmax(h_GP) ≈ 2875.614f0

    # uncertainty
    std_GP = NCDataset(rec_file)["std_uncertainty"][:,:]
    # test that missing values are the same as for surface
    @test findall(ismissing.(std_GP)) == findall(ismissing.(h_GP))
    @test missmin(std_GP) ≈ 0.08829f0 && missmax(std_GP) ≈ 38.07248f0
end


#######################
# GP CROSS-VALIDATION #
#######################

only_atm = true

# get I_no_ocean, (de-)standardization functions and variogram from pre-processing
df_all = CSV.read(csv_preprocessing, DataFrame)
coords_obs = [GrIS1980s_DEM.F.([x,y]) for (x,y) in zip(df_all.x, df_all.y)]
dict   = load(jld2_preprocessing)
@unpack I_no_ocean, idx_aero, gamma, gamma_error, href_file = dict
@unpack standardize, destandardize = GrIS1980s_DEM.get_stddization_fcts(jld2_preprocessing)
varg = GrIS1980s_DEM.get_var(gamma)
varg_error = GrIS1980s_DEM.get_var(gamma_error, nVmax=1)
kernel_signal = varg_to_kernel(varg)
kernel_error  = varg_to_kernel(varg_error)

# make geotable
i_atm = findall(df_all.source .== "atm")
geotable_all = GrIS1980s_DEM.make_geotable(df_all.dh_detrend, df_all.x, df_all.y)
geotable_atm = GrIS1980s_DEM.make_geotable(df_all.dh_detrend[i_atm], df_all.x[i_atm], df_all.y[i_atm])

# add measurement uncertainty for AeroDEM / ATM
GrIS1980s_DEM.add_sigma_obs!(df_all, standardize)

# create sets of training and test data
ℓ         = 2e5
δl        = 1.4e5
geotable_partition = only_atm ? geotable_atm : geotable_all
blocks    = partition(geotable_partition, BlockPartition(ℓ))
ids_test  = indices(blocks)
ids_train = []
for i in eachindex(ids_test)
    if only_atm ids_test[i] = i_atm[ids_test[i]] end
    x_min, x_max = extrema(first.(coords_obs[ids_test[i]]))
    y_min, y_max = extrema(last.(coords_obs[ids_test[i]]))
    i_big = findall(x_min-δl .<= df_all.x .<= x_max+δl .&& y_min-δl .<= df_all.y .<= y_max+δl)
    push!(ids_train, setdiff(i_big, ids_test[i])) # indices that are in i_big but not in ids_test[i]
end

println("GP cross-validation...")
function evaluate_fun(i_train, i_test)
    m_pred, σ_pred = GrIS1980s_DEM.do_GP(coords_obs[i_train], GrIS1980s_DEM.F.(df_all.dh_detrend[i_train]), coords_obs[i_test], kernel_signal, kernel_error, df_all.sigma_obs[i_train], var=true)
    return m_pred, σ_pred
end
difs, sigmas = GrIS1980s_DEM.step_through_folds(ids_train, ids_test, evaluate_fun, df_all.dh_detrend)

# calculate distance to closest observation
dists = nearest_neighb_distance_from_cv(ids_train, ids_test, df_all.x, df_all.y)

# use this distance to bin cross-validation errors
dh_binned, bin_centers = GrIS1980s_DEM.bin_equal_bin_size(vcat(dists...), difs, 7)

@testset "GP cross-validation" begin
    @test !any(isnan.(mean.(dh_binned)))
    @test !any(isnan.(std.(dh_binned)))
    @test maximum(mean.(dh_binned)) ≈ 0.56306f0
    @test maximum(std.(dh_binned)) < 0.9
end


######################
# SVD RECONSTRUCTION  (at 600m resolution with less training data from model realizations)
######################

# working directory back to main dir, to use 600 m input data
cd(pkgdir(GrIS1980s_DEM))
grd_SVD = 600

# define inputs
outline_shp_file = joinpath(pkgdir(GrIS1980s_DEM), "data", "outline", "gris-outline-imbie-1980_updated.shp")
csv_preprocessing, jld2_preprocessing = GrIS1980s_DEM.prepare_obs(grd_SVD, outline_shp_file) # is directed to test/ folder, doesn't recognize the other files
model_realization_files = readdir(joinpath(pkgdir(GrIS1980s_DEM), "data", "model_realizations"), join=true)[1:6]
λ                       = 1e5
r                       = 50

# do reconstruction
rec_file, dict_file = SVD_reconstruction(λ, r, grd_SVD, model_realization_files, csv_preprocessing, jld2_preprocessing)

# test
h_SVD = NCDataset(rec_file)["surface"][:,:]
_, aero_gr_file = GrIS1980s_DEM.create_aerodem(grd_SVD, 150, outline_shp_file)
h_aero = NCDataset(aero_gr_file)["Band1"][:,:]
@testset "SVD reconstruction" begin
    dif = h_aero .- h_SVD
    @test missmax(h_SVD) ≈ 3543.9678f0
    @test mean(abs.(dif[.!ismissing.(dif)])) ≈ 19.5218396
end

# remove output files again
rm(rec_file, recursive=true, force=true)
rm(dict_file, recursive=true, force=true)


########################
# SVD CROSS-VALIDATION #
########################

only_atm  = false

main_output_dir = joinpath("test", "output","validation")
mkpath(main_output_dir)

# get I_no_ocean, (de-)standardization functions and variogram from pre-processing
df_all = CSV.read(csv_preprocessing, DataFrame)
dict   = load(jld2_preprocessing)
@unpack I_no_ocean, idx_aero, href_file = dict

# load datasets, take full SVD (to be truncated with different rs later) for different numbers of training data files
fname = joinpath(main_output_dir, "SVD_components_g$(grd_SVD)_test.jld2")
GrIS1980s_DEM.prepare_model(model_realization_files, I_no_ocean, fname) # read in model data and take svd to derive "eigen ice sheets

# get I_obs and i_atm
@unpack data_mean = load(fname)
x_data, I_obs, i_atm = GrIS1980s_DEM.prepare_obs_SVD(grd_SVD, csv_preprocessing, I_no_ocean, data_mean, main_output_dir)

# create geotable (for partitioning)
x = NCDataset(aero_gr_file)["x"][:]
y = NCDataset(aero_gr_file)["y"][:]
x_Iobs   = x[get_ix.(I_no_ocean[I_obs],length(x))]
y_Iobs   = y[get_iy.(I_no_ocean[I_obs],length(x))]
coords_Iobs = [[ix,iy] for (ix,iy) in zip(x_Iobs,y_Iobs)]
geotable = georef(nothing, coords_Iobs)
geotable_atm = georef(nothing, coords_Iobs[i_atm])

# create sets of training and test data
ℓ    = 2e5 # radius for leaving out Block
logℓ = round(log(10,ℓ),digits=1)
geotable_partition = only_atm ? geotable_atm : geotable
blocks = partition(geotable_partition, BlockPartition(ℓ))
ids_test  = indices(blocks)
ids_train = []
for i in eachindex(ids_test)
    if only_atm ids_test[i] = i_atm[ids_test[i]] end
    id_train    = Vector(1:length(I_obs))
    deleteat!(id_train, ids_test[i])
    push!(ids_train, id_train)
end

# give λ and r values to loop through
λs        = [1e4]
rs        = [5, 20]

function do_validation_and_save(f)
    # load data
    @unpack U, Σ, data_mean, nfiles = load(f)
    UΣ = U*diagm(Σ)
    norms_UΣ = [norm(rw) for rw in eachrow(UΣ)]

    # load datasets, take full SVD (to be truncated with different rs later)
    # x_data, _, _ = GrIS1980s_DEM.prepare_obs_SVD(grd_SVD, csv_preprocessing, I_no_ocean, data_mean, main_output_dir)

    function predict_vals(λ, r, i_train, i_test, x_data, I_obs, UΣ)
        _, x_rec = GrIS1980s_DEM.solve_optim(UΣ, I_obs[i_train], r, λ, x_data[i_train])
        return x_rec[I_obs[i_test]], nothing
    end

    # loop through λ and r values
    m_difs  = [Float32[] for i in eachindex(λs), j in eachindex(rs)]
    for (iλ,λ) in enumerate(λs)
        for (ir,r) in enumerate(rs)
            if r > size(UΣ,2) continue end
            logλ = round(log(10, λ),digits=1)
            println("r = $r, logλ = $logλ")
            evaluate_fun(i_train,i_test) = predict_vals(λ, r, i_train, i_test, x_data, I_obs, UΣ)
            difs, _ = GrIS1980s_DEM.step_through_folds(ids_train, ids_test, evaluate_fun, x_data)
            m_difs[iλ,ir]  = difs
        end
    end
    idx = vcat(ids_test...)
    # save
    to_save = (; norms_UΣ, Σ, dict, grd_SVD, λs, rs, m_difs, idx=I_no_ocean[I_obs[idx]], nfiles, only_atm, method="SVD")
    dest = get_cv_file_SVD(grd_SVD, nfiles; logℓ, only_atm)
    jldsave(dest; to_save...)
end
do_validation_and_save(fname)

dict_file = get_cv_file_SVD(grd_SVD, length(model_realization_files); logℓ, only_atm)
@unpack rs, m_difs = load(dict_file)
@testset "SVD cross-validation" begin
    @test std(m_difs[argmax(rs)]) < std(m_difs[argmin(rs)])
    @test mean(m_difs[2]) ≈ -3.03519f0
    @test std(m_difs[2]) ≈ 26.686956f0
end

rm(dict_file, recursive=true, force=true)


######################
# REMOVE FILES AGAIN #
######################

# remove all test outputs
for file in glob("**/*$(grd).*", joinpath(pkgdir(GrIS1980s_DEM), "data"))
    rm(file, recursive=true, force=true)
end
for file in glob("**/*.*", joinpath(pkgdir(GrIS1980s_DEM), "test", "output"))
    rm(file, recursive=true, force=true)
end
rm(joinpath(pkgdir(GrIS1980s_DEM), "data", "ATM", "atm_dh_interpolated.csv"), recursive=true, force=true)
rm(joinpath(pkgdir(GrIS1980s_DEM), "data", "ATM", "ATM_nadir2seg_all_aligned.csv"), recursive=true, force=true)
rm(joinpath(pkgdir(GrIS1980s_DEM), "data", "ATM", "dh_ref_minus_atm_blocksize_$(blockspacing)m.csv"), recursive=true, force=true)
