using NetCDF, Plots, LinearAlgebra

file = "../Downloads/ex_gris_g1800m_v3a_rcp_26_id_0_2008_2300.nc"

ncinfo(file)

vel = ncread(file, "velsurf_mag") # doesn't seem to conatin any useful information
usurf = ncread(file, "usurf") # surface elevation
ice_mass = ncread(file, "ice_mass")

testgap = ones(size(usurf[:,:,1]))
testgap[350:600,600:1150] .= NaN
testgap[290:520,800:1240] .= NaN
testgap[300:400,200:800]  .= NaN
testgap[400:500,430:600]  .= NaN
u1 = usurf[:,:,1] .* testgap
heatmap(u1')

x1, x2, x3 = [vec(usurf[:,:,t]) for t = [80,90,100,200,250,292] ];
x = [x1 x2 x3];

# optimization problem

U, S, V = svd(x);