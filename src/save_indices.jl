using NetCDF

function get_indices(obs, res)
    # load imbie mask
    imbie_mask = ncread("imbie_mask_g$(res).nc", "Band1")
    no_ocean_mask   = findall((vec(obs) .> 0.0) .|| (vec(imbie_mask) .== 1))
    # get indices where there is data and ice, with respect to ice_mask
    R          = obs[no_ocean_mask]  # vector
    I_obs      = findall(R .> 0.0)
    return no_ocean_mask, I_obs
end
