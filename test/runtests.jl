using svd_IceSheetDEM
using Test, LinearAlgebra, NetCDF, NCDatasets

cd(@__DIR__)  # set working directory to where file is located

const gr = 4000

template_file       = "testdata/testtemplate_g$(gr).nc"
shp_file            = "testdata/testshape.shp"

rm("data/", recursive=true, force=true)
rm("output/", recursive=true, force=true)

###################################################
# testing the data downloading and pre-processing #
###################################################

# bedmachine
bedmachine_file_gr = create_bedmachine_grid(gr, template_file)
bedmachine_path    = splitdir(bedmachine_file_gr)[1]
# aerodem
aero_150_file, aero_gr_file = create_aerodem(;gr, shp_file, bedmachine_path, kw="1981")    # only download the DEMs from 1981 while testing to save some space
# imbie mask
imbie_mask_file = create_imbie_mask(;gr, shp_file, sample_path=aero_150_file)
# atm
atm_file = create_atm_grid(gr, bedmachine_file_gr, "1994.05.21/")

@testset "bedmachine" begin
    ds = NCDataset(bedmachine_file_gr)
    @test all(["mask","geoid","bed","surface","thickness"] .∈ (keys(ds),))
    mask = ds["mask"][:]
    @test sum(mask.==1) .== 31060 && sum(mask.==2) .== 102509 && sum(mask.==3) .== 205 && sum(mask.==4) .== 6894
    @test maximum(ds["geoid"][:]) == 64 && minimum(ds["geoid"][:]) == 0
    @test maximum(ds["bed"][:]) == 3106.2961f0 && minimum(ds["bed"][:]) == -5521.8154f0
    close(ds)
end

@testset "aerodem" begin
    ds150                 = ncread(aero_150_file, "surface")
    dsgr                  = ncread(aero_gr_file, "Band1")
    @test sum(ds150 .> 0) == 2094088 && maximum(ds150)  ≈ 3685.5481
    @test sum(dsgr  .> 0) == 7024    && maximum(dsgr)   ≈ 3252.3076
end

@testset "imbie mask" begin
    imb = ncread(imbie_mask_file, "Band1")
    @test sum(imb .== 1) == 9496
end

@testset "atm" begin
    atm = ncread(atm_file, "surface")
    @test sum(atm .> 0) == 633
    @test maximum(atm) == 2385.68f0
end

###############################
# testing the problem solving #
###############################

F           = Float32
λ           = 1e5
r           = 10^3
rec_file    = solve_lsqfit(F, λ, r, gr, imbie_mask_file, bedmachine_file_gr, [template_file], aero_gr_file)
rec_bm_file = create_reconstructed_bedmachine(rec_file, bedmachine_file_gr)

@testset "solve least square fit" begin
    rec         = ncread(rec_file, "surface")
    @test maximum(rec) ≈ 3086.8433f0 && rec[202:205, 79] ≈ Float32[492.5, 440.9, 530.0, 557.6] && sum(rec .> 0) == 5247
    bm          = NCDataset(rec_bm_file)
    @test all(["mask","bed","surface","thickness","polar_stereographic","x","y"] .∈ (keys(bm),))
    mask = bm["mask"][:]
    @test sum(mask.==1) == 31060 && sum(mask.==2) == 4999 && sum(mask.==3) == 212
    @test isapprox(maximum(bm["surface"][:]), 3106, atol=1)
    close(bm)
end

####################################
# testing the uncertainty analysis #
####################################
