library(raster)
library(xts)

#
PISCOts <- seq(as.Date("1981-01-01"),
               as.Date("2016-12-31"), by = "day")
PISCOts <- xts::xts(1:length(PISCOts), PISCOts)

# 
tmax_normal <- brick(
  lapply(1:12, function(x){
    raster(file.path("data", "processed", "PISCOt_DS",
            paste("tmax_",
                  ifelse(x < 10, formatC(x, 1, flag = "0"), as.character(x)),
                  ".nc", sep = "")))
    })
  )

tmin_normal <- brick(
  lapply(1:12, function(x){
    raster(file.path("data", "processed", "PISCOt_DS",
                     paste("tmin_",
                           ifelse(x < 10, formatC(x, 1, flag = "0"), as.character(x)),
                           ".nc", sep = "")))
  })
)

parallel::mclapply(1:length(PISCOts),
                   function(i){
          
  anom_tmax <- raster(file.path(".", "data", "processed", "PISCOt", "anomalies",
                                paste("tmax_", time(PISCOts[i]), ".nc", sep = "")))

  anom_tmax <- disaggregate(anom_tmax, 10, method = "bilinear")
  #anom_tmax <- mask(anom_tmax, tmax_normal[[1]])
  anom_tmax <- resample(anom_tmax, tmax_normal[[1]])

  PISCOtmax <- tmax_normal[[ as.numeric(format(time(PISCOts[i]), "%m")) ]] - anom_tmax

  writeRaster(PISCOtmax, file.path(".", "data", "processed", "PISCOt_DS", "values", "tmax",
                                   paste("tmax_", time(PISCOts[i]), ".nc", sep = "")),
              overwrite = TRUE)

  ##
  
  anom_tmin <- raster(file.path(".", "data", "processed", "PISCOt", "anomalies",
                                paste("tmin_", time(PISCOts[i]), ".nc", sep = "")))

  anom_tmin <- disaggregate(anom_tmin, 10, method = "bilinear")
  #anom_tmin <- mask(anom_tmin, tmin_normal[[1]])
  anom_tmin <- resample(anom_tmin, tmin_normal[[1]])

  PISCOtmin <- tmin_normal[[ as.numeric(format(time(PISCOts[i]), "%m")) ]] - anom_tmin

  writeRaster(PISCOtmin, file.path(".", "data", "processed", "PISCOt_DS", "values", "tmin",
                                   paste("tmin_", time(PISCOts[i]), ".nc", sep = "")),
              overwrite = TRUE)

  }, mc.cores = 8)

#fun::shutdown()