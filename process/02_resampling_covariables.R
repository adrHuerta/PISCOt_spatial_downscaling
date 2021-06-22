library(raster)

#
PISCOtmax_normal <- raster("./data/processed/PISCOt/PISCOtmax_normal.nc")

#
lst_day <- brick("./data/raw/LST/LST_DAY.nc")
lst_night <- brick("./data/raw/LST/LST_NIGHT.nc")
z <- brick("./data/raw/Z/DEM.nc")[[1]]

#
lst_day_res <- resample(lst_day, PISCOtmax_normal[[1]])
lst_night_res <- resample(lst_night, PISCOtmax_normal[[1]])
z_res <- resample(z, PISCOtmax_normal[[1]])

# same resolution and dims?
extent(PISCOtmax_normal) == extent(lst_day_res)
extent(PISCOtmax_normal) == extent(lst_night_res)
extent(PISCOtmax_normal) == extent(z_res)

dim(PISCOtmax_normal) == dim(lst_day_res[[1]])
dim(PISCOtmax_normal) == dim(lst_night_res[[1]])
dim(PISCOtmax_normal) == dim(z_res)

#
writeRaster(lst_day_res, "./data/processed/LST/LST_DAYres.nc", overwrite = TRUE,
            datatype = 'FLT4S', force_v4 = TRUE, compression = 7)

writeRaster(lst_night_res, "./data/processed/LST/LST_NIGHTres.nc", overwrite = TRUE,
            datatype = 'FLT4S', force_v4 = TRUE, compression = 7)

writeRaster(z_res, "./data/processed/Z/Zres.nc", overwrite = TRUE,
            datatype = 'FLT4S', force_v4 = TRUE, compression = 7)

