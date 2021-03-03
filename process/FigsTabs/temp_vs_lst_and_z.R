library(raster)
"%>%" = magrittr::`%>%`

#
tmax <- brick("./data/processed/PISCOt/PISCOtmax_normal.nc")
tmin <- brick("./data/processed/PISCOt/PISCOtmin_normal.nc")

#
lst_day <- brick("./data/processed/LST/LST_DAYres.nc")
lst_night <- brick("./data/processed/LST/LST_NIGHTres.nc")
z <- raster("./data/processed/Z/Zres.nc")

#
tmax_pairs <- brick(mean(tmax), mean(lst_day), z)
names(tmax_pairs) <- c("Tmax", "LST_día", "Z")

jpeg("./output/FigsTabs/temp_vs_lst_z_tmax.jpg", 
     width = 950, height = 950, res = 200)
pairs(tmax_pairs, cex = .1)
dev.off()

jpeg("./output/FigsTabs/temp_vs_lst_z_tmax_ag10.jpg", 
     width = 950, height = 950, res = 200)
pairs(aggregate(tmax_pairs, 10), cex = .1)
dev.off()

jpeg("./output/FigsTabs/temp_vs_lst_z_tmax_de10.jpg", 
     width = 950, height = 950, res = 200)
pairs(disaggregate(tmax_pairs, 10, method = "bilinear"), cex = .1)
dev.off()
#

tmin_pairs <- brick(mean(tmin), mean(lst_night), z)
names(tmin_pairs) <- c("Tmin", "LST_noche", "Z")

jpeg("./output/FigsTabs/temp_vs_lst_z_tmin.jpg", 
     width = 950, height = 950, res = 200)
pairs(tmin_pairs, cex = .1)
dev.off()





