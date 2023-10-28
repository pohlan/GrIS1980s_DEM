using svd_IceSheetDEM, Glob


# different method for testing that it's the same result

# rec = ncread("output/bedmachine1980_reconstructed_g600.nc", "surface")
# x   = ncread("output/bedmachine1980_reconstructed_g600.nc", "x")
# y   = ncread("output/bedmachine1980_reconstructed_g600.nc", "y")
# dx  = x[2]-x[1]
# dy  = y[2]-y[1]
# geoid = ncread("data/bedmachine/bedmachine_g600.nc", "geoid")

# sc = CSV.read("data/ATM/ATM_nadir2seg_all.csv", DataFrame)
# coords = [[sc[i,3], sc[i,2]] for i in eachindex(sc[:,2])]
# coords_proj = AG.reproject(coords, ProjString("+proj=longlat +datum=WGS84 +no_defs"), EPSG(3413))
# x_proj = first.(coords_proj)
# y_proj = last.(coords_proj)

# npts = zeros(Int, size(rec))
# atm_grid = zeros(Float32, size(rec))
# z_val = zeros(Float32, size(rec))
# idxs = findall(rec .!= 0.0)
# @showprogress for ci in idxs[1:10000]
#     ix, iy = ci.I
#     pts_in_reach = findall((abs.(x[ix] .- x_proj) .< 0.5dx) .&& (abs.(y[iy] .- y_proj) .< 0.5dy))
#     npts[ix,iy] = length(pts_in_reach)
#     if npts[ix,iy] >= 3
#         rs = sqrt.((x_proj[pts_in_reach].-x[ix]).^2 .+ (y_proj[pts_in_reach].-y[iy]).^2)
#         z_val[ix,iy] = sum(sc[pts_in_reach,4] ./ rs)  /  sum(1 ./ rs) - geoid[ix,iy]
#     end
# end




















# dhdt_files = String[]
# for (root, dirs, files) in walkdir("data/ATM/dh_dt/")
#     for file in files
#         if startswith(file, "IDHDT4") && endswith(file, ".csv")
#             println(joinpath(root, file)) # path to files
#             push!(dhdt_files, joinpath(root, file))
#         end
#     end
# end
# csv_mosaic_to_one_netcdf(atm_path, dhdt_files, "dhdt_all")

# r((lat1, lon1), (lat2, lon2)) = 6.378e3^2*(2 - 2*cos(lat1)*cos(lat2)*cos(lon1-lon2) - 2*sin(lat1)*sin(lat2))

# function plot_it()
#     Pl.figure()
#     for endyr in 1998:2018   #["1998", "2002", "2003", "2007", "2009"]
#         for startyr in 1993:1994
#             fi = findall(occursin.("$endyr-$startyr", dhdt_files))
#             if !isempty(fi)
#                 dhdt = CSV.read(dhdt_files[fi], DataFrame)
#                 Pl.scatter(dhdt[:,2],dhdt[:,1],5,dhdt[:,4],cmap="bwr", clim=(-1,1))
#                 println("$endyr-$startyr")
#             end
#         end
#     end
#     Pl.colorbar()
#     Pl.savefig("scatter_all.jpg")
# end
# plot_it()
