using TSVD, NetCDF, Statistics, LinearAlgebra, Glob, PyPlot, Printf

include("read_in.jl")
model_files = glob("data/usurf_ex_gris_g1800m_v5_RAGIS_id_*.nc") # 60 time steps in total, 27 files
Data, x, y, t, nx, ny = training_data(files=5:5,tsteps=1:60;model_files)

# total volume and change
V_tot = [sum(Data[i,:]) for i in 1:size(Data,1)]
dV_tot = diff(V_tot)

# compute SVD
q = 10
U, S, V = tsvd(Data, q)

# compute volume change for each eigenglacier
dVs = zeros(q,1)
for i = 1:q
    eigs   = U[:,i] * S[i] * V[:,i]'
    V_eig  = [sum(eigs[i,:]) for i in 1:size(Data,1)]
    dVs[i] = sum(diff(V_eig))
end

# contribution of eigenglacier vol. change to total vol. change
r_ind =               dVs *100 ./ sum(dV_tot)
r_cum = cumsum(dVs,dims=1)*100 ./ sum(dV_tot)

@printf(" %1.1f%% of the volume change is captured in the first eigenglacier.\n", r_cum[1])
@printf(" %1.1f%% of the volume change is captured in the first + second eigenglaciers.\n", r_cum[2])
