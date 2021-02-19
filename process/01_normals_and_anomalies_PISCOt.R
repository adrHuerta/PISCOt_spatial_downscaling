library(raster)
library(xts)

#a
PISCOtmax <- brick("./data/raw/PISCOt/PISCOdtx_v1.1.nc")
PISCOtmin <- brick("./data/raw/PISCOt/PISCOdtn_v1.1.nc")

#
PISCOts <- xts::xts(1:nlayers(PISCOtmax), seq(as.Date("1981-01-01"),
                                              as.Date("2016-12-31"), by = "day"))

# normals 1981-2010
PISCOts_1981_2010 <- PISCOts["1981/2010"]

PISCOtmax_normal <- brick(
  lapply(1:12, function(z){
    PISCOidx <- PISCOts_1981_2010[.indexmon(PISCOts_1981_2010) %in% (z-1)]
    mean(PISCOtmax[[PISCOidx]])
    })
  )

PISCOtmax_normal <- setNames(PISCOtmax_normal,1:12)

PISCOtmin_normal <- brick(
  lapply(1:12, function(z){
    PISCOidx <- PISCOts_1981_2010[.indexmon(PISCOts_1981_2010) %in% (z-1)]
    mean(PISCOtmin[[PISCOidx]])
  })
)

PISCOtmin_normal <- setNames(PISCOtmin_normal,1:12)

writeRaster(PISCOtmax_normal, "./data/processed/PISCOt/PISCOtmax_normal.nc")
writeRaster(PISCOtmin_normal, "./data/processed/PISCOt/PISCOtmin_normal.nc")

# anomalies
for(i in 1:length(seq_along(PISCOts))){
  
  ni <- as.numeric(PISCOts[i])
  mi <- (.indexmon(PISCOts[i]) + 1)
  tmax_normal <- PISCOtmax_normal[[paste("X", mi, sep = "")]]
  tmin_normal <- PISCOtmin_normal[[paste("X", mi, sep = "")]]
  
  anom_tmax <- tmax_normal - PISCOtmax[[ni]]
  anom_tmin <- tmin_normal - PISCOtmin[[ni]]
  
  writeRaster(anom_tmax, file.path(".", "data", "processed", "PISCOt", "anomalies",
                              paste("tmax_", time(PISCOts[i]), ".nc", sep = "")))
  writeRaster(anom_tmin, file.path(".", "data", "processed", "PISCOt", "anomalies",
                              paste("tmin_", time(PISCOts[i]), ".nc", sep = "")))
}
