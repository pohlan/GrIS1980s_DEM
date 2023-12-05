using svd_IceSheetDEM
using Test, LinearAlgebra, NetCDF, NCDatasets

cd(@__DIR__)  # set working directory to where file is located

#####################################
# testing the randomized svd (rsvd) #
#####################################

@testset "test rsvd with m = $m, n = $n" for
    (m, n) in ((100, 100),
               (1000, 150),
               (50, 900))
    A = randn(m,n)
    U, S, V = svd(A)
    # for k in [2]
        Ur, Sr, Vr = rsvd(A, min(m,n))
        @test norm(Sr - S) < sqrt(eps(eltype(A)))
        @test norm(abs.(transpose(Ur)*U) - I) < sqrt(eps(eltype(A)))
        @test norm(abs.(transpose(Vr)*V) - I) < sqrt(eps(eltype(A)))
    # end
end
@testset "truncated rsvd" begin
    A = [1. 2 3; 4 5 6; 7 8 9]
    U, S, V = svd(A)
    for k = 2:3
        Ur, Sr, Vr = rsvd(A, k)
        @test norm(Sr - S[1:k]) < sqrt(eps(eltype(A)))
        @test norm(abs.(transpose(Ur)*U[:,1:k]) - I) < sqrt(eps(eltype(A)))
        @test norm(abs.(transpose(Vr)*V[:,1:k]) - I) < sqrt(eps(eltype(A)))
    end
end

#######################################
# testing preprocessing gdal routines #
#######################################

const gr = 4000

imbie_path          = "data/gris-imbie-1980/"
aerodem_path        = "data/aerodem/"
aero_150            = aerodem_path*"aerodem_rm-filtered_geoid-corr_g150.nc"
aero_gr             = aerodem_path*"aerodem_rm-filtered_geoid-corr_g$(gr).nc"
bedmachine_path     = "data/bedmachine/"
bedmachine_file_gr  = bedmachine_path*"bedmachine_g$(gr).nc"
template_file       = "testdata/testtemplate_g$(gr).nc"
shp_file            = "testdata/testshape.shp"
imbie_mask          = imbie_path * "imbie_mask_g$(gr).nc"

@testset "bedmachine" begin
    rm(bedmachine_path, recursive=true, force=true)
    create_bedmachine_grid(gr, bedmachine_path, template_file)
    ds = NCDataset(bedmachine_file_gr)
    @test all(["mask","geoid","bed","surface","thickness"] .∈ (keys(ds),))
    mask = ds["mask"][:]
    @test sum(mask.==1) .== 144181 && sum(mask.==2) .== 111489 && sum(mask.==3) .== 536 && sum(mask.==4) .== 6450
    @test maximum(ds["geoid"][:]) == 64 && minimum(ds["geoid"][:]) == 0
    @test maximum(ds["bed"][:]) == 3106.2961f0 && minimum(ds["bed"][:]) == -5521.8154f0
    close(ds)
end

@testset "aerodem" begin
    rm(aerodem_path, recursive=true, force=true)
    create_aerodem(aerodem_path, shp_file, bedmachine_path, "1981")    # only download the DEMs from 1981 while testing to save some space
    gdalwarp(aero_150; gr, srcnodata="0.0", dest=aero_gr)
    ds150 = ncread(aero_150, "surface")
    dsgr  = ncread(aero_gr , "Band1")
    @test sum(ds150 .> 0) == 2094088 && maximum(ds150) == 3685.548095703125
    @test sum(dsgr  .> 0) == 7024    && maximum(dsgr)  == 3252.3076520647323
end

@testset "imbie mask" begin
    rm(imbie_path, recursive=true, force=true)
    create_imbie_mask(gr; imbie_path, imbie_shp_file=shp_file, sample_path=aero_150)
    imb = ncread(imbie_mask, "Band1")
    @test sum(imb .== 1) == 9496
end

###############################
# testing the problem solving #
###############################

@testset "solve least square fit" begin
    rm("../output/", recursive=true, force=true)
    F           = Float32
    λ           = 1e5
    r           = 1e3
    rec_file    = solve_lsqfit(F, λ, r, gr, imbie_mask, [template_file], aero_gr)
    rec         = ncread(rec_file, "surface")
    @test maximum(rec) ≈ 3042.707f0 && rec[200:203, 80] ≈ Float32[62.8, 62.8, 89.6, 147.5] && sum(rec .> 0) == 9144
    create_reconstructed_bedmachine(rec_file, bedmachine_file_gr)
    bm          = NCDataset("output/bedmachine1980_reconstructed_g$(gr).nc")
    @test all(["mask","bed","surface","thickness","polar_stereographic","x","y"] .∈ (keys(bm),))
    mask = bm["mask"][:]
    @test sum(mask.==1) == 144046 && sum(mask.==2) == 5824 && sum(mask.==3) == 135
    @test isapprox(maximum(bm["surface"][:]), 3043, atol=1)
    close(bm)
end
