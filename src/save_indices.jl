using NetCDF

function get_indices(obs, res)
    # load bedmachine mask
    mask            = ncread("data/pism_Greenland_" * res * "m_mcb_jpl_v2023_RAGIS_ctrl.nc", "mask")
    nx, ny          = size(mask)
    not_ocean_mask  = mask.!= 0. .&& mask .!= 4  # mask for model grid point that are neither ocean nor outside Greenland

    # extend the not_ocean_mask a little bit so that it does not cut off glaciers that extended further out in the 1980's
    not_ocean_extended = copy(not_ocean_mask)
    l = 10
    for ix = l+1:size(not_ocean_extended,1)-l
        for iy = l+1:size(not_ocean_extended,2)-l
            if any(not_ocean_mask[ix-l:ix+l,iy-l:iy+l])
                not_ocean_extended[ix,iy] = true
            end
        end
    end

    I_not_ocean     = findall(reshape(not_ocean_extended,nx*ny))                    # non-ocean indices
    # mask_ice_nested = reshape(mask,nx*ny)[I_not_ocean] .== 2
    no_bedrock_nested = reshape(mask,nx*ny)[I_not_ocean] .!= 1

    # get indices where there is data and ice, with respect to non-ocean-mask
    R         = reshape(obs, length(obs))[I_not_ocean]
    I_marg    = findall(R .!= -9999 .&& R .!= 0.0 .&& no_bedrock_nested)
    # get indices of the interior where there is ice but no data, with respect to non-ocean-mask
    bool_intr = (obs.==-9999 .|| obs .== 0.0)
    for iy = 1:ny
        ix = 1
        while (obs[ix,iy] == -9999 || obs[ix,iy] == 0.0)
            bool_intr[ix,iy] = false
            ix += 1
            if ix > nx
                break
            end
        end
        ix = size(obs,1)
        while (obs[ix,iy] == -9999 || obs[ix,iy] == 0.0)
            bool_intr[ix,iy] = false
            ix -= 1
            if ix < 1
                break
            end
        end
    end
    I_intr = findall(reshape(bool_intr,nx*ny,1)[I_not_ocean] .&& no_bedrock_nested) ####### mask = ice, data = no interior of ice sheet

    return I_not_ocean, I_marg, I_intr
end
