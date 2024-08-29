using JLD2, Plots, LaTeXStrings

fig_dir = "output/SVD_reconstruction/figures"


########################
# plot eigenmodes      #
########################

# define how many modes to plot
n_modes = 4

# load
d = load("output/SVD_reconstruction/SVD_components.jld2")
U = d["U"][:,1:n_modes]
I_no_ocean = d["I_no_ocean"]

nx, ny = 2640, 4560   # hard-coded for 600m

b = zeros(nx*ny, n_modes)
b[I_no_ocean, :] .= U
b[b .== 0] .= NaN

Plots.scalefontsizes()
Plots.scalefontsizes(3.0)

colors = cgrad(:RdBu)

# plot eigenmodes
ps = [heatmap(reshape(b[:,i], nx,ny)', clims=(-2e-3, 2e-3), cmap=:RdBu, showaxis=false, top_margin = -20Plots.mm, leftmargin=-30Plots.mm, grid=false, title= latexstring("u_$i"), cbar=false, aspect_ration=true) for i in 1:n_modes]
p_panels = plot(ps..., size=(2100,900), aspect_ratio=1, left_margin=[0Plots.mm 0Plots.cm], right_margin=[0Plots.mm 0Plots.mm], bottom_margin=-9Plots.mm, layout=Plots.grid(1,length(ps)))

# add a colorbar separately (inspired by https://discourse.julialang.org/t/set-custom-colorbar-tick-labels/69806/4)
xx = range(0,1,500)
zz = zero(xx)' .+ xx
p_c = heatmap(xx, xx, zz, ticks=false, ratio=25, legend=false, fc=colors, lims=(0,1), clims=(0,1), framestyle=:box, right_margin=20Plots.mm, top_margin=10Plots.mm, bottom_margin=15Plots.mm)
yticks  = [0.05, 0.5, 0.95]
ticktxt = ["negative", "0", "positive"]
[annotate!(1.5, yi, text(ti, 20, "Computer Modern", :left)) for (yi,ti) in zip(yticks,ticktxt)]

# plot again everything together
layout = @layout [a{0.9w,1.0h} b{0.7h}]
plot(p_panels, p_c; layout, bottom_margin=-40Plots.mm, size=(2100,600), top_margin=10Plots.mm)
Plots.savefig(joinpath(fig_dir, "eigenmodes.png"))



########################
# plot v vectors       #
########################

# plot(d["v_rec"], label="")
# plot!(d["V"][20,:], label="")


