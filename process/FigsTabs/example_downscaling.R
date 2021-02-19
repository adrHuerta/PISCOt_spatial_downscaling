library(raster)
"%>%" = magrittr::`%>%`

#
tmax <- brick("./data/processed/PISCOt/PISCOtmax_normal.nc")
lst_day <- brick("./data/processed/LST/LST_DAYres.nc")
z <- brick("./data/processed/Z/Zres.nc")[[1]]

#tmax
# at 5
i = 2
brick_model <- aggregate(brick(tmax[[i]], lst_day[[i]], z), 3)
names(brick_model) <- c("Temp", "LST", "Z")

data_model <- rasterToPoints(brick_model, fun=NULL, spatial=TRUE)
data_model <- data_model[complete.cases(data_model@data),]

bw <- GWmodel::bw.gwr(Temp ~ LST + Z, data = data_model,
                      approach = "AICc",
                      kernel = "bisquare",
                      adaptive = TRUE,
                      longlat=TRUE)

localstats1_lst <- GWmodel::gwss(data_model, vars = c("Temp", "LST"), bw = 24)
localstats1_z <- GWmodel::gwss(data_model, vars = c("Temp", "Z"), bw = 24)

brick(rasterFromXYZ(localstats1_lst$SDF[, "Corr_Temp.LST"]),
      rasterFromXYZ(localstats1_z$SDF[, "Corr_Temp.Z"])) %>%
        setNames(c("cor(Tmax,LST)", "cor(Tmax,Z)")) 

rasterFromXYZ(localstats1_lst$SDF[, "Corr_Temp.LST"]) %>%
        spplot()

rasterFromXYZ(localstats1_z$SDF[, "Corr_Temp.Z"]) %>%
        spplot()

brick(,
      rasterFromXYZ(localstats1_z$SDF[, "Corr_Temp.Z"])) %>%
        setNames(c("cor(Tmax,LST)", "cor(Tmax,Z)"))



model_sel <- GWmodel::model.selection.gwr(DeVar= c("Temp"), 
                             InDeVars = c("LST","Z"),
                             approach = "AICc",
                             data = data_model,
                             bw = bw,
                             kernel = "bisquare",
                             adaptive = TRUE,
                             longlat=TRUE) 
sorted_models <- GWmodel::gwr.model.sort(model_sel, numVars = 2,
                                         ruler.vector = model_sel[[2]][,2])
model_list <- sorted_models[[1]]                             
GWmodel::model.view.gwr(c("Temp"), c("LST","Z"), model_list)

bw_sel <- c(1393, 869, 543, 344, 218, 143, 94, 66, 46, 36, 27, 24, 20, 24)
AICc_sel <- c(6460.027, 6197.105, 5752.631, 5053.354,
              4348.042, 3715.732, 3220.133, 2902.301, 
              2648.609, 2522.297, 2425.73, 2408.637, 
              2630.694, 2408.637) 


jpeg("./output/FigsTabs/model_selection.jpg", 
     width = 1000, height = 700, res = 200)
plot(sorted_models[[2]][, 2], type = "o", col = "black",
     xlab = "Modelo", ylab = "AICc", xaxt = "n", 
     xlim = c(0.75, 3.25),
     ylim = c(1000, 4500))
axis(1, at=c(1, 2, 3), labels=c("Temp ~ LST", "Temp ~ Z", "Temp ~ LST+Z"))
dev.off()

jpeg("./output/FigsTabs/bw_selection.jpg", 
     width = 1000, height = 700, res = 200)
plot(bw_sel, AICc_sel, type = "o",
     xlab = "Bandwith", ylab = "AICc", col = "black")
dev.off()


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


Intercept <- rasterFromXYZ(gwr_model$SDF[, "Intercept"])
LST_coef <- rasterFromXYZ(gwr_model$SDF[, "LST"])
Z_coef <- rasterFromXYZ(gwr_model$SDF[, "Z"])

GWmodel_coef <- brick(Intercept, LST_coef, Z_coef)
crs(GWmodel_coef) <- crs(brick_model)
GWestimate <- GWmodel_coef[[1]] + (GWmodel_coef[[2]]*brick_model[[2]]) + (GWmodel_coef[[3]]*brick_model[[3]])
GW_residual <- brick_model[[1]] - GWestimate

GW_output1 <- brick(Intercept, LST_coef, Z_coef, GW_residual)
names(GW_output1) <- c("Intercept", "LST_coef", "Z_coef", "GW_residual")

##
LST_day <- brick("./data/raw/LST/LST_DAY.nc")
LST_night <- brick("./data/raw/LST/LST_NIGHT.nc")
Z <- brick("./data/raw/Z/DEM.nc")

r0 <- LST_day[[1]]
r1 <- LST_night[[1]]
r2 <- Z[[1]]
r0 <- crop(r0, r2)
r2 <- crop(r2, r0)
r1 <- crop(r1, r2)

r <- brick(r0, r1, r2)[[1]]
#

GW_output_res <- disaggregate(GW_output1, 10, method = "bilinear")
GW_output_res <- crop(GW_output_res, r)
GW_output_res <- resample(GW_output_res, r)

Intercpt <- spplot(GW_output_res[[1]], at = seq(-50, 50, 5),  col.regions=terrain.colors(500),
                   main = "Intercepto")
LST_coef <- spplot(GW_output_res[[2]], at = seq(-1.5, 1.5, .1),  col.regions=terrain.colors(500),
                   main = "LST")
Z_coef <- spplot(GW_output_res[[3]]*100, at = seq(-0.1, .1, .005)*100,  col.regions=terrain.colors(500),
                 main = "Z")
Resdl_model <- spplot(GW_output_res[[4]], at = seq(-2, 2, .01),  col.regions=terrain.colors(500),
                      main = "Residual")
library(latticeExtra)
require(gridExtra)
jpeg("./output/FigsTabs/model_coef_int_lst.jpg", 
     width = 1000, height = 800, res = 200)
grid.arrange(Intercpt,LST_coef, ncol=2)
dev.off()

jpeg("./output/FigsTabs/model_coef_z_residual.jpg", 
     width = 1000, height = 800, res = 200)
grid.arrange(Z_coef, Resdl_model, ncol=2)
dev.off()


# example 14-02-2009

PISCOtmax <- brick("./data/raw/PISCOt/PISCOdtx_v1.1.nc")

#
PISCOts <- xts::xts(1:nlayers(PISCOtmax), seq(as.Date("1981-01-01"),
                                              as.Date("2016-12-31"), by = "day"))

PISCOtmax_DS <- brick("./data/processed/PISCOt_DS/values/tmax/tmax_2000-02-14.nc") 

p3 <- spplot(crop(PISCOtmax[[PISCOts["2000-02-14"]]],  extent(-80, -75, -10, -5)), at = seq(-5, 35, 1)) 
p4 <- spplot(crop(PISCOtmax_DS[[1]],  extent(-80, -75, -10, -5)), at = seq(-5, 35, 1))
p1 <- spplot(crop(PISCOtmax[[PISCOts["2000-02-14"]]],  extent(-75, -70, -15, -10)), at = seq(-5, 35, 1)) 
p2 <- spplot(crop(PISCOtmax_DS[[1]],  extent(-75, -70, -15, -10)), at = seq(-5, 35, 1))

c(p3, p4, p1, p2)
     
adr <- crop(PISCOtmax[[PISCOts["2000-02-14"]]],  extent(-75, -70, -15, -10))
adr2 <- resample(adr, crop(PISCOtmax_DS[[1]],  extent(-75, -70, -15, -10)), method = "bilinear")

brick(adr2, PISCOtmax_DS[[1]])
spplot((PISCOtmax_DS[[1]]-adr2))

p1<-rasterVis::levelplot(crop(PISCOtmax[[PISCOts["2000-02-14"]]],  extent(-75, -70, -15, -10)), contour=TRUE,
                     at = seq(-4, 36, 4), margin = list(draw=F),
                     scales=list(draw=FALSE), xlab = "", ylab = "")
p2<-rasterVis::levelplot(crop(PISCOtmax_DS[[1]],  extent(-75, -70, -15, -10)), contour=TRUE,
                     at = seq(-4, 36, 4), margin = list(draw=F),
                     scales=list(draw=FALSE), xlab = "", ylab = "")

p3<-rasterVis::levelplot(crop(PISCOtmax[[PISCOts["2000-02-14"]]],  extent(-80, -75, -10, -5)), contour=TRUE,
                         at = seq(-4, 36, 4), margin = list(draw=F),
                         scales=list(draw=FALSE), xlab = "", ylab = "")
p4<-rasterVis::levelplot(crop(PISCOtmax_DS[[1]],  extent(-80, -75, -10, -5)), contour=TRUE,
                         at = seq(-4, 36, 4), margin = list(draw=F),
                         scales=list(draw=FALSE), xlab = "", ylab = "")

png("./output/FigsTabs/comparison_example.png", 
     width = 1200, height = 1200, res = 200)
print(c(p3, p4, p1, p2))
dev.off()

