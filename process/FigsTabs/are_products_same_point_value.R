library(zoo)
library(raster)

qc_data <- readRDS("./data/processed/PISCOt_v1-2/QC_GF_HG_data.RDS")
xyz <- qc_data$xyz 

# PISCOt v1.1
PISCOt <- brick("./data/raw/PISCOt/PISCOdtn_v1.1.nc")
grid_cells_1 <- extract(PISCOt[[1]],
                        xyz[c(100, 50), c("LON", "LAT")],
                        cellnumbers = TRUE)[, "cells"]

ts_example_1 <- t(PISCOt[grid_cells_1])

# PISCO v1.1 downscaled
grid_PISCOt_scale_reduced_tmin <- raster("/home/waldo/Documentos/Repos/netcdf_files/PISCOt_SD/tmin/tmin_1981-01.nc")
grid_PISCOt_scale_reduced_tmax <- raster("/home/waldo/Documentos/Repos/netcdf_files/PISCOt_SD/tmax/tmax_1981-01.nc")

grid_cells_2 <- extract(grid_PISCOt_scale_reduced_tmin,
                        xyz[c(100, 50), c("LON", "LAT")],
                        cellnumbers = TRUE)[, "cells"]


PISCOt_v1_1_ds_tmin <- dir("/home/waldo/Documentos/Repos/netcdf_files/PISCOt_SD/tmin",
    full.names = TRUE)

PISCOt_v1_1_ds_tmax <- dir("/home/waldo/Documentos/Repos/netcdf_files/PISCOt_SD/tmax",
                           full.names = TRUE)

ts_example_2 <- do.call("rbind",
                        parallel::mclapply(PISCOt_v1_1_ds_tmin,
                                           function(x){
                                             t(brick(x)[grid_cells_2])
                                             }, mc.cores = 6))

jpeg("./output/FigsTabs/comparison_01.jpg", 
     width = 1000, height = 600, res = 150)
matplot(cbind(ts_example_2[,2], ts_example_1[,2]), cex = 0.25, type = "p", pch = 19,
        ylab = "", main = "80.65058 -3.939097  142")
dev.off()

jpeg("./output/FigsTabs/comparison_02.jpg", 
     width = 1000, height = 600, res = 150)
matplot(cbind(ts_example_2[,1], ts_example_1[,1]), cex = 0.25, type = "p", pch = 19,
        ylab = "", main = "-78.80512 -6.379639 2668")
dev.off()
## 

grid_cells_tmin <- extract(grid_PISCOt_scale_reduced_tmin,
                        xyz[, c("LON", "LAT")],
                        cellnumbers = TRUE)[, "cells"]

grid_cells_tmax <- extract(grid_PISCOt_scale_reduced_tmax,
                           xyz[, c("LON", "LAT")],
                           cellnumbers = TRUE)[, "cells"]

ts_tmin <- do.call("rbind",
                        parallel::mclapply(PISCOt_v1_1_ds_tmin,
                                           function(x){
                                             t(brick(x)[grid_cells_tmin])
                                           }, mc.cores = 6))

ts_tmax <- do.call("rbind",
                   parallel::mclapply(PISCOt_v1_1_ds_tmax,
                                      function(x){
                                        t(brick(x)[grid_cells_tmax])
                                      }, mc.cores = 6))

dim(qc_data$values$tmax["1981/2016"])
dim(ts_tmax)

qc_data$xyz$tmax_r <- NA
qc_data$xyz$tmin_r <- NA
qc_data$xyz$tmax_b <- NA
qc_data$xyz$tmin_b <- NA

for(i in 1:ncol(ts_tmax)){
  qc_data$xyz$tmax_r[i] <- cor(qc_data$values$tmax["1981/2016"][, i],
                               ts_tmax[, i], use = "pairwise.complete.obs")
  
  qc_data$xyz$tmin_r[i] <- cor(qc_data$values$tmin["1981/2016"][, i],
                               ts_tmin[, i], use = "pairwise.complete.obs")
  
  qc_data$xyz$tmax_b[i] <- hydroGOF::pbias(obs = as.numeric(qc_data$values$tmax["1981/2016"][, i]),
                               sim = as.numeric(ts_tmax[, i]), na.rm = TRUE)
  
  qc_data$xyz$tmin_b[i] <- hydroGOF::pbias(obs = as.numeric(qc_data$values$tmin["1981/2016"][, i]),
                               sim = as.numeric(ts_tmin[, i]), na.rm = TRUE)

}

qc_data_exp <- qc_data$xyz[complete.cases(qc_data$xyz), ]
qc_data_exp <- transform(qc_data_exp,
                         tmax_r = cut(tmax_r, 
                                      breaks = c(-Inf, .3, .4, .5, .6, .7, .8, .9, Inf),
                                      labels = c("< .3",  ".3 - .4", ".4 - .5", ".5 - .6", ".6 - .7", ".7 - .8", ".8 - .9","> .9"),
                                      right = FALSE))

qc_data_exp <- transform(qc_data_exp,
                         tmin_r = cut(tmin_r, 
                                      breaks = c(-Inf, .3, .4, .5, .6, .7, .8, .9, Inf),
                                      labels = c("< .3",  ".3 - .4", ".4 - .5", ".5 - .6", ".6 - .7", ".7 - .8", ".8 - .9","> .9"),
                                      right = FALSE))

qc_data_exp <- transform(qc_data_exp,
                         tmax_b = cut(tmax_b, 
                                      breaks = c(-Inf, -20, -10, -5, 5, 10, 20, Inf),
                                      labels = c("< -20",  "-20 - -10", "-10 - -5", "-5 - 5", "5 - 10", "10 - 20","> 20"),
                                      right = FALSE))

qc_data_exp <- transform(qc_data_exp,
                         tmin_b = cut(tmin_b, 
                                      breaks = c(-Inf, -20, -10, -5, 5, 10, 20, Inf),
                                      labels = c("< -20",  "-20 - -10", "-10 - -5", "-5 - 5", "5 - 10", "10 - 20","> 20"),
                                      right = FALSE))

library(lattice)
library(sp)

shp_peru = file.path(".", "data", "raw", "vectorial", "Departamentos.shp") %>% 
  raster::shapefile() 

sp::SpatialPointsDataFrame(coords = qc_data_exp[, c("LON", "LAT")],
                           data = qc_data_exp[, c(9,10, 11, 12)],
                           proj4string = CRS("+init=epsg:28992")) -> sp_data

jpeg("./output/FigsTabs/comparison_tmax_r.jpg", 
     width = 1000, height = 1000, res = 150)
sp::spplot(sp_data[, "tmax_r"]) + 
  latticeExtra::layer(sp.polygons(shp_peru, fill = NA, col = "gray50"), under = TRUE, superpose = FALSE) 
dev.off()

jpeg("./output/FigsTabs/comparison_tmin_r.jpg", 
     width = 1000, height = 1000, res = 150)
sp::spplot(sp_data[, "tmin_r"]) + 
  latticeExtra::layer(sp.polygons(shp_peru, fill = NA, col = "gray50"), under = TRUE, superpose = FALSE)
dev.off()


jpeg("./output/FigsTabs/comparison_tmin_b.jpg", 
     width = 1000, height = 1000, res = 150)
sp::spplot(sp_data[, "tmin_b"]) + 
  latticeExtra::layer(sp.polygons(shp_peru, fill = NA, col = "gray50"), under = TRUE, superpose = FALSE)
dev.off()

jpeg("./output/FigsTabs/comparison_tmax_b.jpg", 
     width = 1000, height = 1000, res = 150)
sp::spplot(sp_data[, "tmax_b"]) + 
  latticeExtra::layer(sp.polygons(shp_peru, fill = NA, col = "gray50"), under = TRUE, superpose = FALSE)
dev.off()