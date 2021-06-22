library(raster)
library(gstat)

# grid of 0.01
LST_day <- brick("./data/raw/LST/LST_DAY.nc")
LST_night <- brick("./data/raw/LST/LST_NIGHT.nc")
Z <- brick("./data/raw/Z/DEM.nc")

r0 <- LST_day[[1]]
r1 <- LST_night[[1]]
r2 <- Z[[1]]
r0 <- mask(r0, r2)
r2 <- mask(r2, r0)
r1 <- mask(r1, r2)

r <- brick(r0, r1, r2)[[1]]

LST_day <- mask(LST_day, r)
LST_night <- mask(LST_night, r)
Z <- mask(Z[[1]], r)

for(i in 1:12){
  
  GW_file <- file.path("data", "processed", "GW_model",
                       paste("tmax_", "GWoutput_",
                             ifelse(i < 10, formatC(i, 1, flag = "0"), as.character(i)),
                             ".nc", sep = ""))
  
  GW_output <- brick(GW_file)
  names(GW_output) <- c("Intercept", "LST_coef", "Z_coef", "GW_residual")
  
  GW_output_res <- disaggregate(GW_output, 10, method = "bilinear")
  #GW_output_res <- mask(GW_output_res, r)
  GW_output_res <- resample(GW_output_res, r)
  
  GW_modelClim <- GW_output_res[[1]] + GW_output_res[[2]]*LST_day[[1]] + GW_output_res[[3]]*Z
  GW_modelClim <- GW_output_res[[4]] + GW_modelClim
  
  writeRaster(GW_modelClim, 
              file.path("data", "processed", "PISCOt_DS",
                        paste("tmax_",
                              ifelse(i < 10, formatC(i, 1, flag = "0"), as.character(i)),
                              ".nc", sep = "")),
              overwrite = TRUE,
              datatype = 'FLT4S', force_v4 = TRUE, compression = 7)
}


for(i in 1:12){
  
  GW_file <- file.path("data", "processed", "GW_model",
                       paste("tmin_", "GWoutput_",
                             ifelse(i < 10, formatC(i, 1, flag = "0"), as.character(i)),
                             ".nc", sep = ""))
  
  GW_output <- brick(GW_file)
  names(GW_output) <- c("Intercept", "LST_coef", "Z_coef", "GW_residual")
  
  GW_output_res <- disaggregate(GW_output, 10, method = "bilinear")
  #GW_output_res <- mask(GW_output_res, r)
  GW_output_res <- resample(GW_output_res, r)
  
  GW_modelClim <- GW_output_res[[1]] + GW_output_res[[2]]*LST_night[[1]] + GW_output_res[[3]]*Z
  GW_modelClim <- GW_output_res[[4]] + GW_modelClim
  
  writeRaster(GW_modelClim, 
              file.path("data", "processed", "PISCOt_DS",
                        paste("tmin_",
                              ifelse(i < 10, formatC(i, 1, flag = "0"), as.character(i)),
                              ".nc", sep = "")),
              overwrite = TRUE,
              datatype = 'FLT4S', force_v4 = TRUE, compression = 7)
}
