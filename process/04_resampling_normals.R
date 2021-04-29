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
              overwrite = TRUE)
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
              overwrite = TRUE)
}




# r0 <- brick("./data/raw/LST/LST_DAY.nc")[[1]]
# r1 <- crop(brick("./data/raw/Z/DEM.nc")[[1]], r0)
# r0 <- crop(r0, r1)
# r <- brick(r0, r1)[[1]]
# crs(r) <- crs(brick("./data/processed/GW_model/tmax_GWoutput_01.nc"))
# 
# # projection(r) <- projection(tmin)
# r2 <- aggregate(r, 7)
# 
# 
# #
# GW_exp <- brick("./data/processed/GW_model/tmax_GWoutput_01.nc")
# names(GW_exp) <- c("Intercept", "LST_coef", "Z_coef", "GW_residual")
# 
# 
# 
# 
# 
# 
# 
# 
# 
# GW_exp <- brick("./data/processed/GW_model/tmax_GWoutput_01.nc")
# names(GW_exp) <- c("Intercept", "LST_coef", "Z_coef", "GW_residual")
# GW_exp <- rasterToPoints(GW_exp$Intercept, fun=NULL, spatial=TRUE)
# GW_exp <- 
# 
# gs <- gstat::gstat(formula = Intercept ~ 1, 
#                    locations = GW_exp[complete.cases(GW_exp@data$Intercept), ],
#                    set = list(idp=2))
# system.time(idw <- interpolate(r2, gs))
# plot(idw)












# 
# library(fields)
# m <- fastTps(coordinates(GW_exp), GW_exp$Intercept, theta =3)
# tps <- interpolate(r, m)
# 
# 
# gs <- gstat::gstat(formula=Intercept~1, locations=GW_exp)
# idw <- interpolate(r, gs)
# 
# plot(idw)
# plot(LST_coef)
# g <- as(r2, 'SpatialPixelsDataFrame')
# 
# m <- automap::autofitVariogram(formula=LST~1, input_data=gwr_model$SDF[, "LST"],
#                                fix.values = c(0,NA,NA))
# k <- gstat::gstat(formula=LST~1, locations=gwr_model$SDF[, "LST"], model=m$var_model)
# kp <- raster::predict(k, g)
# 
# library(parallel)
# # Calculate the number of cores
# no_cores <- 3
# 
# # Initiate cluster (after loading all the necessary object to R environment: meuse, meuse.grid, m)
# cl <- makeCluster(no_cores)
# 
# parts <- split(x = 1:length(g), f = 1:no_cores)
# 
# clusterExport(cl = cl, varlist = c("gwr_model", "g", "m", "parts"), envir = .GlobalEnv)
# clusterEvalQ(cl = cl, expr = c(library('sp'), library('gstat'), library("raster")))
# 
# parallelX <- parLapply(cl = cl, X = 1:no_cores, 
#                        fun = function(x){
#                          k <- gstat::gstat(formula=LST~1, locations=gwr_model$SDF[, "LST"], model=m$var_model)
#                          raster::predict(k, g[parts[[x]], ])
#                        })
