############################
# Plot SVD modes, Figure 3 #
############################

# define which modes to plot
i_modes = [1,2,3,100]
# load
d = load(joinpath("output", "reconstructions", "SVD_components_g$(grd)_rec.jld2"))
U = d["U"][:,i_modes]
# load into array of full map
nx, ny = length(x), length(y)
b = zeros(nx*ny, length(i_modes))
b[I_no_ocean, :] .= U
b[b .== 0] .= NaN
# plotting attributes
scalefontsizes()
scalefontsizes(3.0)
colors = cgrad(:RdBu)
# plot modes
ps = [heatmap(reshape(b[:,i], nx,ny)', clims=(-2e-3, 2e-3), cmap=:RdBu, showaxis=false, top_margin = -20Plots.mm, leftmargin=-30Plots.mm, grid=false, title= latexstring("u_{$(i_modes[i])}"), cbar=false, aspect_ration=true) for i in 1:length(i_modes)]
p_panels = plot(ps..., size=(2100,900), aspect_ratio=1, left_margin=[0Plots.mm 0Plots.cm], right_margin=[0Plots.mm 0Plots.mm], bottom_margin=-9Plots.mm, layout=Plots.grid(1,length(ps)))
# add a colorbar separately (inspired by https://discourse.julialang.org/t/set-custom-colorbar-tick-labels/69806/4)
xx = range(0,1,500)
zz = zero(xx)' .+ xx
p_c = heatmap(xx, xx, zz, ticks=false, ratio=25, legend=false, fc=colors, lims=(0,1), clims=(0,1), framestyle=:box, right_margin=20Plots.mm, top_margin=10Plots.mm, bottom_margin=15Plots.mm)
yticks  = [0.05, 0.5, 0.95]
ticktxt = ["negative", "0", "positive"]
[annotate!(1.5, yi, text(ti, 20, "Computer Modern", :left)) for (yi,ti) in zip(yticks,ticktxt)]
# plot everything together
layout = @layout [a{0.9w,1.0h} b{0.7h}]
plot(p_panels, p_c; layout, bottom_margin=-40Plots.mm, size=(2100,600), top_margin=10Plots.mm)
savefig(joinpath(fig_dir_main, "Figure3.png"))


########################
# Plot v vectors       #
########################

fig_dir_others = joinpath("output", "reconstructions", "figures")
mkpath(fig_dir_others)

v_rec = d["v_rec"]
bar(d["v_rec"][1:100], label="", linecolor=:cornflowerblue, color=:cornflowerblue, size=attr.size, margin=attr.margin, xlabel=L"Mode index $i$", ylabel=L"$\mathbf{v}_\mathrm{rec,\,i}$", grid=false)
savefig(joinpath(fig_dir_others, "v_rec.png"))
