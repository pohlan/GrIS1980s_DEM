using svd_IceSheetDEM
import ArchGDAL as AG

# grid
gr = 1200

# load datasets
aeroDEM    = shortread("data/aerodem_raw/aerodem_rm-filtered_geoid-corr_g$(gr).nc")
surfaceDEM = shortread("output/rec_lambda_1e5_g$(gr)_r377.nc")
bedDEM     = shortread("data/bedm_bed_g$(gr).nc")
bedm_mask  = shortread("data/bedm_mask_g$(gr).nc")
ice_mask   = (surfaceDEM .> 0.0) .&& (surfaceDEM .> bedDEM)

# calculate floating mask
ρw            = 1030
ρi            = 917
Pw            = - ρw * bedDEM
Pi            = ρi * (surfaceDEM - bedDEM)
floating_mask = (Pw .> Pi) .&& (ice_mask)  # floating where water pressure > ice pressure at the bed

# calculate mask
new_mask = zeros(Int, size(bedm_mask))
new_mask[ice_mask] .= 2
new_mask[bedm_mask .== 1] .= 1 # bedrock
new_mask[floating_mask]   .= 3

# calculate ice thickness
h_ice = zeros(size(surfaceDEM))
h_ice[ice_mask] .= surfaceDEM[ice_mask] - bedDEM[ice_mask]
h_ice[floating_mask] .= surfaceDEM[floating_mask] ./  (1-ρi/ρw)

# save to netcdf file
template_path = "data/aerodem_raw/aerodem_rm-filtered_geoid-corr_g1200.nc"
dest = "output/bedmachine1980_reconstructed_test.nc"
template_dataset = AG.read(template_path)
    AG.create(
        dest,
        driver = AG.getdriver(template_dataset),
        width  = AG.width(template_dataset),
        height = AG.height(template_dataset),
        nbands = 4,
        dtype  = Float32,
        options = ["-co", "COMPRESS=DEFLATE", "-co", "ZLEVEL=6"] # reduces the file size
    ) do raster
        AG.write!(raster, surfaceDEM, 1)
        AG.write!(raster, bedDEM, 2)
        AG.write!(raster, h_ice, 3)
        AG.write!(raster, new_mask, 4)
        AG.setgeotransform!(raster, AG.getgeotransform(template_dataset))
        AG.setproj!(raster, AG.getproj(template_dataset))
    end

# couldn't figure out how to name a layer with ArchGDAL, this hack probably needs a fix...
for (n, varb) in enumerate(["surface", "bed", "thickness", "mask"])
    run(`ncrename -v $("Band"*string(n)),$varb $dest`)
end
