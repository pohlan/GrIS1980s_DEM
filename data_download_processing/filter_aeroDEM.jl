using NetCDF, Statistics, ImageFiltering
import ArchGDAL as AG

# read in data
bedm_orig   = ncread("../data/BedMachineGreenland-v5.nc", "surface")
bedm_rot    = bedm_orig[:,end:-1:1]                                # has to be flipped to match aerodem
aerodem     = ncread("../data/aerodem_g150m_geoid_corrected_1978_1987_mean.nc", "surface")

# find indices where aerodem has observations
i_aero      = aerodem .!= -9999.0

# take difference
dif         = NaN .* ones(size(aerodem))
dif[i_aero] = bedm_rot[i_aero] .- aerodem[i_aero]

# apply median filter
flt_median = x -> sum(isnan.(x))/length(x) .> 0.8 ? NaN : x[ceil(Int,end/2),ceil(Int,end/2)]
flt_slop   = x -> any(abs.(diff(x,dims=1)) .> 30) || any(abs.(diff(x,dims=2)) .> 30) ? NaN : x[ceil(Int,end/2),ceil(Int,end/2)]
dif_fmed   = mapwindow(flt_median, dif, (5,5))
dif_clean  = mapwindow(flt_slop, dif_fmed, (3,3))


# create aerodem with data only where dif_clean has data
i_clean    = .!isnan.(dif_clean)
aerodem_clean = zeros(size(aerodem))
aerodem_clean[i_clean] .= aerodem[i_clean]

# create netcdf file by copying attributes from original file
sample_dataset = AG.read("../data/aerodem_g150m_geoid_corrected_1978_1987_mean.nc")
AG.create(
        "../data/aerodem_filtered_g150m_geoid.nc",
        driver = AG.getdriver(sample_dataset),
        width  = AG.width(sample_dataset),
        height = AG.height(sample_dataset),
        nbands = 1,
        dtype  = Float32
    ) do raster
        AG.write!(raster, aerodem_clean[:,end:-1:1], 1)
        AG.setgeotransform!(raster, AG.getgeotransform(sample_dataset))
        AG.setproj!(raster, AG.getproj(sample_dataset))
    end
