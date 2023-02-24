using JLD2, NetCDF

mkpath("data/")


resolutions = ["1200", "1800"]
for res in resolutions

    # load bedmachine mask
    mask            = ncread("data/pism_Greenland_" * res * "m_mcb_jpl_v2023_RAGIS_ctrl.nc", "mask")
    nx, ny          = size(mask)
    not_ocean_mask  = reshape(mask.!= 0. .&& mask .!= 4, nx*ny)  # mask for model grid point that are neither ocean nor outside Greenland
    I_not_ocean     = findall(not_ocean_mask)                    # non-ocean indices
    mask_ice_nested = reshape(mask,nx*ny)[I_not_ocean] .== 2     

    # load the aerodem observations
    obs_file ="data/aerodem_g" * res * "m_geoid_corrected_1978_1987_mean.nc"
    obs      = ncread(obs_file, "surface_altitude")
    R        = reshape(obs, length(obs))[I_not_ocean]

    # get indices where there is data and ice, with respect to non-ocean-mask
    I_marg    = findall(R .!= -9999 .&& mask_ice_nested)
    # get indices of the interior where there is ice and but no data, with respect to non-ocean-mask
    bool_intr = obs.==-9999
    for iy = 1:ny
        ix = 1
        while obs[ix,iy] .== -9999
            bool_intr[ix,iy] = false
            ix += 1
            if ix > nx
                break
            end
        end
        ix = size(obs,1)
        while obs[ix,iy] .== -9999
            bool_intr[ix,iy] = false
            ix -= 1
            if ix < 1
                break
            end
        end
    end
    I_intr = findall(reshape(bool_intr,nx*ny,1)[I_not_ocean] .&& mask_ice_nested) ####### mask = ice, data = no interior of ice sheet

    # save
    jldsave("data/indices_aerodem_g" * res * "m.jld2"; I_not_ocean, I_marg, I_intr)
end
