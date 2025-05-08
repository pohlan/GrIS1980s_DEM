using GrIS1980s_DEM
using Test, LinearAlgebra, NetCDF, NCDatasets, CSV, DataFrames

cd(@__DIR__)  # set working directory to where file is located

const grd = 4000

test_svd_training_file      = joinpath("testdata", "testtemplate_g$(grd).nc")
outline_shp_file            = joinpath("testdata", "testshape.shp")

rm("data", recursive=true, force=true)
rm("output", recursive=true, force=true)

###################################################
# testing the data downloading and pre-processing #
###################################################

# bedmachine
bedmachine_original, bedmachine_file_gr = GrIS1980s_DEM.create_bedmachine_grid(grd)
# GrIMP v2
tiles = ["2_0", "2_1", "3_0", "3_1", "3_2", "4_1", "4_2", "5_2"]
ref_coreg_file_ellips, ref_coreg_file_geoid, ref_file   = GrIS1980s_DEM.create_grimpv2(grd, 150, bedmachine_original, kw=tiles)
# aerodem
aero_150_file, aero_gr_file = GrIS1980s_DEM.create_aerodem(grd, 150, outline_shp_file, bedmachine_original, ref_coreg_file_geoid, kw="1981")    # only download the DEMs from 1981 while testing to save some space
# imbie mask
outline_mask_file = GrIS1980s_DEM.create_outline_mask(grd, outline_shp_file, aero_150_file)
# atm
blockspacing=grd
lines = ["93.06.23", "93.06.24", "93.07.01", "93.07.09"]
atm_dh_file = GrIS1980s_DEM.get_atm_dh_file(ref_coreg_file_ellips, ref_coreg_file_geoid, outline_shp_file, blockspacing, kw=lines)

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
    @test missmax(ds) ≈ 3184.565f0
    @test sum(.!ismissing.(ds)) == 35396
    # check that it's not "upside down"
    ny_half = fld(size(ds,2),2)
    @test sum(.!ismissing.(ds[:,(ny_half+1):end])) == 0    # The data should all be in the southern half
end

@testset "aerodem" begin
    ds150                 = NCDataset(aero_150_file)["surface"][:,:]
    ds150                 = nomissing(ds150,0.0)
    dsgr                  = NCDataset(aero_gr_file )["Band1"][:,:]
    @test missmax(ds150)  ≈ 3637.929f0
    @test sum(.!ismissing.(dsgr)) == 2971 && missmax(dsgr) ≈ 3096.081f0
    # check that it's not "upside down"
    ny_half = fld(size(dsgr,2),2)
    sum(.!ismissing.(dsgr[:,(ny_half+1):end])) == 0    # The data should all be in the southern half
end

@testset "outline mask" begin
    outline_mask = NCDataset(outline_mask_file)["Band1"][:,:]
    @test missum(outline_mask .== 1) == sum(.!ismissing.(outline_mask)) == 9050
    # check that it's not "upside down"
    ny_half = fld(size(outline_mask,2),2)
    sum(.!ismissing.(outline_mask[:,(ny_half+1):end])) == 0    # The data should all be in the southern half
end

@testset "atm" begin
    df_atm = CSV.read(atm_dh_file, DataFrame)
    @test size(df_atm) == (2555,4)
    @test all(extrema(df_atm.x) .≈ (-303464f0, 638898f0))
    @test all(extrema(df_atm.y) .≈ (-3.294195f6, -1.814327f6))
    @test sum(.!ismissing.(df_atm.dh)) == 1498  # the ones outside of the outline are missing because there is no h_ref
    @test missmax(df_atm.h_ref) ≈ 3176f0
    @test missmin(df_atm.dh) ≈ -1580.9f0 && missmax(df_atm.dh) ≈ 49.8f0
end
