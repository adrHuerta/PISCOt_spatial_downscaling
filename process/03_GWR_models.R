library(GWmodel)
library(raster)

lst_day <- brick("./data/processed/LST/LST_DAYres.nc")
lst_night <- brick("./data/processed/LST/LST_NIGHTres.nc")
z <- brick("./data/processed/Z/Zres.nc")[[1]]
tmax <- brick("./data/processed/PISCOt/PISCOtmax_normal.nc")
tmin <- brick("./data/processed/PISCOt/PISCOtmin_normal.nc")

for(i in 1:12){
  
  #tmax
  # at 5
  brick_model <- aggregate(brick(tmax[[i]], lst_day[[i]], z), 3)
  names(brick_model) <- c("Temp", "LST", "Z")
  
  data_model <- rasterToPoints(brick_model, fun=NULL, spatial=TRUE)
  data_model <- data_model[complete.cases(data_model@data),]
  
  bw <- GWmodel::bw.gwr(Temp ~ LST + Z, data = data_model,
                        approach = "AICc",
                        kernel = "bisquare",
                        adaptive = TRUE,
                        longlat=TRUE)
  
  # dmat_dis <- GWmodel::gw.dist(dp.locat= as.matrix(as.data.frame(data_model)[, c(4,5)]), longlat=TRUE,
  #                              rp.locat = as.matrix(as.data.frame(data_model)[, c(4,5)]))
  # 
  # GWmodel::gwr.mixed(Temp ~ LST + Z, 
  #                    data = data_model, 
  #                    bw = bw,  
  #                    kernel = "bisquare", 
  #                    adaptive = TRUE,
  #                    longlat=TRUE,
  #                    regression.points = newpts,
  #                    dMat = dmat_dis)
  # 
  # original
  brick_model <- brick(tmax[[i]],lst_day[[i]], z)
  names(brick_model) <- c("Temp", "LST", "Z")
  
  data_model <- rasterToPoints(brick_model, 
                               fun=NULL, spatial=TRUE)
  data_model <- data_model[complete.cases(data_model@data),]
  newpts <- rasterToPoints(brick_model[[1]]); newpts <- newpts[, -3]
  
  gwr_model <- tryCatch(GWmodel::gwr.basic(Temp ~ LST + Z, 
                                           data = data_model, 
                                           bw = bw,  
                                           kernel = "bisquare", 
                                           adaptive = TRUE,
                                           longlat=TRUE,
                                           regression.points = newpts, 
                                           parallel.method = "cluster"),
                        error = function(x){ NA })
  
  if(all(is.na(gwr_model))){
    
    brick_model <- aggregate(brick(tmax[[i]], lst_day[[i]], z), 2)
    names(brick_model) <- c("Temp", "LST", "Z")
    
    data_model <- rasterToPoints(brick_model, 
                                 fun=NULL, spatial=TRUE)
    data_model <- data_model[complete.cases(data_model@data),]
    newpts <- rasterToPoints(brick_model[[1]]); newpts <- newpts[, -3]
    
    gwr_model <- tryCatch(GWmodel::gwr.basic(Temp ~ LST + Z, 
                                             data = data_model, 
                                             bw = bw,  
                                             kernel = "bisquare", 
                                             adaptive = TRUE,
                                             longlat=TRUE,
                                             regression.points = newpts, 
                                             parallel.method = "cluster"),
                          error = function(x){ NA })
  }
  
  if(all(is.na(gwr_model))){
    
    brick_model <- aggregate(brick(tmax[[i]], lst_day[[i]], z), 3)
    names(brick_model) <- c("Temp", "LST", "Z")
    
    data_model <- rasterToPoints(brick_model, 
                                 fun=NULL, spatial=TRUE)
    data_model <- data_model[complete.cases(data_model@data),]
    newpts <- rasterToPoints(brick_model[[1]]); newpts <- newpts[, -3]
    
    gwr_model <- tryCatch(GWmodel::gwr.basic(Temp ~ LST + Z, 
                                             data = data_model, 
                                             bw = bw,  
                                             kernel = "bisquare", 
                                             adaptive = TRUE,
                                             longlat=TRUE,
                                             regression.points = newpts, 
                                             parallel.method = "cluster"),
                          error = function(x){ NA })
  }
  
  Intercept <- rasterFromXYZ(gwr_model$SDF[, "Intercept"])
  LST_coef <- rasterFromXYZ(gwr_model$SDF[, "LST"])
  Z_coef <- rasterFromXYZ(gwr_model$SDF[, "Z"])
  
  GWmodel_coef <- brick(Intercept, LST_coef, Z_coef)
  crs(GWmodel_coef) <- crs(brick_model)
  GWestimate <- GWmodel_coef[[1]] + (GWmodel_coef[[2]]*brick_model[[2]]) + (GWmodel_coef[[3]]*brick_model[[3]])
  GW_residual <- brick_model[[1]] - GWestimate
  
  GW_output1 <- brick(Intercept, LST_coef, Z_coef, GW_residual)
  #names(GW_output1) <- c("Intercept", "LST_coef", "Z_coef", "GW_residual")
  crs(GW_output1) <- crs(brick_model)
  
  writeRaster(GW_output1, file.path(".", "data", "processed", "GW_model",
                                    paste("tmax_GWoutput_", ifelse(i < 10, formatC(i, 1, flag = "0"), as.character(i)), ".nc", sep = "")), 
              overwrite = TRUE,
              datatype = 'FLT4S', force_v4 = TRUE, compression = 7)
}

for(i in 1:12){
  
  #tmin
  # at 5
  brick_model <- aggregate(brick(tmin[[i]], lst_night[[i]], z), 3)
  names(brick_model) <- c("Temp", "LST", "Z")

  data_model <- rasterToPoints(brick_model, fun=NULL, spatial=TRUE)
  data_model <- data_model[complete.cases(data_model@data),]
  #pairs(data_model@data[, c("Temp", "LST", "Z")])
  bw <- GWmodel::bw.gwr(Temp ~ LST + Z, data = data_model,
                        approach = "AICc",
                        kernel = "bisquare",
                        adaptive = TRUE,
                        longlat=TRUE)
  
  
  # original
  brick_model <- brick(tmin[[i]], lst_night[[i]], z)
  names(brick_model) <- c("Temp", "LST", "Z")
  
  data_model <- rasterToPoints(brick_model, 
                               fun=NULL, spatial=TRUE)
  data_model <- data_model[complete.cases(data_model@data),]
  newpts <- rasterToPoints(brick_model[[1]]); newpts <- newpts[, -3]
  
  gwr_model <- tryCatch(GWmodel::gwr.basic(Temp ~ LST + Z, 
                                  data = data_model, 
                                  bw = bw,  
                                  kernel = "bisquare", 
                                  adaptive = TRUE,
                                  longlat=TRUE,
                                  regression.points = newpts, 
                                  parallel.method = "cluster"),
           error = function(x){ NA })
  
  if(all(is.na(gwr_model))){
    
    brick_model <- aggregate(brick(tmin[[i]], lst_night[[i]], z), 2)
    names(brick_model) <- c("Temp", "LST", "Z")
    
    data_model <- rasterToPoints(brick_model, 
                                 fun=NULL, spatial=TRUE)
    data_model <- data_model[complete.cases(data_model@data),]
    newpts <- rasterToPoints(brick_model[[1]]); newpts <- newpts[, -3]
    
    gwr_model <- tryCatch(GWmodel::gwr.basic(Temp ~ LST + Z, 
                                             data = data_model, 
                                             bw = bw,  
                                             kernel = "bisquare", 
                                             adaptive = TRUE,
                                             longlat=TRUE,
                                             regression.points = newpts, 
                                             parallel.method = "cluster"),
                          error = function(x){ NA })
    }
  
  if(all(is.na(gwr_model))){
    
    brick_model <- aggregate(brick(tmin[[i]], lst_night[[i]], z), 3)
    names(brick_model) <- c("Temp", "LST", "Z")
    
    data_model <- rasterToPoints(brick_model, 
                                 fun=NULL, spatial=TRUE)
    data_model <- data_model[complete.cases(data_model@data),]
    newpts <- rasterToPoints(brick_model[[1]]); newpts <- newpts[, -3]
    
    gwr_model <- tryCatch(GWmodel::gwr.basic(Temp ~ LST + Z, 
                                             data = data_model, 
                                             bw = bw,  
                                             kernel = "bisquare", 
                                             adaptive = TRUE,
                                             longlat=TRUE,
                                             regression.points = newpts, 
                                             parallel.method = "cluster"),
                          error = function(x){ NA })
  }
  
  Intercept <- rasterFromXYZ(gwr_model$SDF[, "Intercept"])
  LST_coef <- rasterFromXYZ(gwr_model$SDF[, "LST"])
  Z_coef <- rasterFromXYZ(gwr_model$SDF[, "Z"])
  
  GWmodel_coef <- brick(Intercept, LST_coef, Z_coef)
  crs(GWmodel_coef) <- crs(brick_model)
  GWestimate <- GWmodel_coef[[1]] + (GWmodel_coef[[2]]*brick_model[[2]]) + (GWmodel_coef[[3]]*brick_model[[3]])
  GW_residual <- brick_model[[1]] - GWestimate
  
  GW_output1 <- brick(Intercept, LST_coef, Z_coef, GW_residual)
  #names(GW_output1) <- c("Intercept", "LST_coef", "Z_coef", "GW_residual")
  crs(GW_output1) <- crs(brick_model)
  
  writeRaster(GW_output1, file.path(".", "data", "processed", "GW_model",
                                    paste("tmin_GWoutput_", ifelse(i < 10, formatC(i, 1, flag = "0"), as.character(i)), ".nc", sep = "")), overwrite = TRUE)
}
