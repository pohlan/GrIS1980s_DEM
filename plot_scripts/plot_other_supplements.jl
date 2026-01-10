######################################
# overlap AeroDEM and ATM, Figure S1 #
######################################

csv_aero_atm = "data/ATM/aero_minus_ATM.csv"
df_dh = CSV.read(csv_aero_atm, DataFrame)
i_nonan = findall(.!ismissing.(df_dh.dh))

scalefontsizes()
scalefontsizes(1.9)
histogram(df_dh.dh[i_nonan], xlims=(-100,100), normalize=:pdf, yticks=false, grid=false,
          fillcolor=:cornflowerblue, linecolor=:cornflowerblue,
          ylabel="Frequency", xlabel="\n"*L"$h_\mathrm{AeroDEM}-h_\mathrm{ATM}\;\;\mathrm{(m)}$", label="",
          size=(900,700), margin=10mm)
savefig(joinpath(fig_dir_main, "fS01.png"))

# map
# scalefontsizes()
# scalefontsizes(2.7)
# plot(multipolygon, fillcolor=:grey, fillalpha=0.2, linewidth=0, label="AeroDEM\ncoverage")
# scatter!(df_dh.x[i_nonan].*1e-3, df_dh.y[i_nonan].*1e-3, marker_z=df_dh.dh[i_nonan], legend_foreground_color=nothing, legend=(0.01,0.4), markersize=4, markerstrokewidth=0, cmap=cgrad(:vik, rev=true), clims=(-50,50), aspect_ratio=1, xaxis=false, yaxis=false, grid=false, size=(1500,1900), label="", margin=10mm, xlims=(-7e5,8e5).*1e-3, ylims=(-3.32e6, -0.78e6).*1e-3)
# savefig(joinpath(fig_dir_main, "fSXX_aero_minus_atm.png"))
