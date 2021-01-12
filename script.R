## ----read_raster------------------------------------------------------
library(raster)
rar <- raster("data/bathy.tif")  
class(rar)


## ----plot_raster------------------------------------------------------
plot(rar)


## ----read_stars-------------------------------------------------------
library(stars)
ras <- read_stars("data/bathy.tif")   # argument `driver` to specify the driver
class(ras)


## ----plot_stars, cache = TRUE-----------------------------------------
plot(ras)


## ----raster_drivers---------------------------------------------------
ra_dr <- rgdal::gdalDrivers()
head(ra_dr)


## ----raster_drivers2--------------------------------------------------
ra_dr[which(ra_dr$name == "GTiff"), ]


## ----stars_drivers----------------------------------------------------
st_dr <- sf::st_drivers(what = "raster")
head(st_dr)


## ----stars_drivers2---------------------------------------------------
st_dr[which(st_dr$name == "GTiff"), ]


## ----create_dir-------------------------------------------------------
dir.create("output", showWarnings = FALSE)


## ----write_raster0----------------------------------------------------
ra_dr$name[ra_dr$create]


## ----write_raster1, cache = TRUE--------------------------------------
writeRaster(rar, filename = "output/rar.gpkg", format = "GPKG", overwrite = TRUE)


## ----write_stars0-----------------------------------------------------
sort(st_dr$name[st_dr$write])


## ----write_stars1, cache = TRUE---------------------------------------
write_stars(ras, dsn = "output/ras.gpkg", driver = "GPKG")




## ----raster_obj-------------------------------------------------------
class(rar)
rar


## ----raster_obj_val---------------------------------------------------
values(rar)


## ----raster_obj_fun---------------------------------------------------
dim(rar)
ncell(rar)
projection(rar)
bbox(rar)
extent(rar)


## ----stars_obj--------------------------------------------------------
class(ras)
ras


## ----stars_val--------------------------------------------------------
class(ras[[1]])
ras[[1]]


## ----stars_obj_fun1---------------------------------------------------
dim(ras)
ncell(ras)
st_bbox(ras)


## ----stars_obj_fun2---------------------------------------------------
ras_crs <- st_crs(ras)
class(ras_crs)
ras_crs


## ----stars_obj_fun3---------------------------------------------------
ras_crs$input
ras_crs$proj4string
ras_crs$epsg


## ----convert1---------------------------------------------------------
rar_c <- st_as_stars(rar) 
class(rar_c)


## ----convert2---------------------------------------------------------
ras_c <- as(ras, "Raster") 
class(ras_c)


## ----bc_rar1----------------------------------------------------------
class(rar[,])
rar[1, 1]
rar[1:10, 11:20]
mean(rar[,]) # or mean(as.matrix(rar))
quantile(rar[,])


## ----bc_rar1_2--------------------------------------------------------
rar2 <- rar # create a copy 
rar2[rar > 0] <- NA  # filter
plot(rar2)


## ----bc_ras1----------------------------------------------------------
# extract value
class(ras[[1]])
ras[[1]][1, 1]
ras[[1]][1:10, 11:20]


## ----bc_ras1_2--------------------------------------------------------
mean(ras[[1]]) 
quantile(ras[[1]])


## ----bc_rar1_3--------------------------------------------------------
ras2 <- ras # create a copy
ras2[[1]][ras[[1]] < units::as_units(0, "m")] <- NA # filter
plot(ras2)




## ----rar_new_crs------------------------------------------------------
rar_crs <- CRS("+proj=merc +a=6378137 +b=6378137 +lat_ts=0.0 +lon_0=0.0 +x_0=0.0 +y_0=0 +k=1.0 +units=m +nadgrids=@null +wktext  +no_defs")


## ----projrar, cache = TRUE--------------------------------------------
rar_t <- projectRaster(rar, crs = rar_crs)
projection(rar_t)
plot(rar)


## ----projrar2, cache = TRUE-------------------------------------------
ras_t <- st_transform(ras, crs = 3857)
st_crs(ras_t)


## ----crop0------------------------------------------------------------
plot(rar)
rect(-65, 45, -60, 50)


## ----crop_rar, cache = TRUE-------------------------------------------
rar_c <- crop(rar, extent(-65, -60, 45, 50))
plot(rar_c)


## ----crop_ras, cache = TRUE-------------------------------------------
ext <- st_bbox(c(xmin = -65, xmax = -60, ymin = 45, ymax = 50), crs = 4326)
ras_c <- st_crop(ras, ext)
plot(ras_c)


## ----stl--------------------------------------------------------------
(stl <- sf::st_read("data/st_laurence.geojson"))


## ----plot_stl, cache = TRUE-------------------------------------------
plot(sf::st_geometry(stl))


## ----rar_m, cache = TRUE----------------------------------------------
rar_m <- mask(rar, stl)
plot(rar_m)


## ----rar_m2, cache = TRUE---------------------------------------------
rar_mi <- mask(rar, stl, inverse = TRUE)
plot(rar_mi)


## ----ras_m, cache = TRUE----------------------------------------------
ras_m <- ras[stl]
plot(ras_m)


## ----ras_m1, cache = TRUE---------------------------------------------
stl_i <- sf::st_difference(sf::st_as_sfc(st_bbox(ras)), stl)
plot(stl_i, col = 2)


## ----ras_m2, cache = TRUE---------------------------------------------
ras_mi <- ras[stl_i]
plot(ras_mi)


## ---- echo = FALSE----------------------------------------------------
countdown::countdown(minutes = 10, seconds = 0)





## ----rar_warp, cache = TRUE-------------------------------------------
# create a raster with the given resolution
mr <- matrix(runif(21*21), 21, 21)
rar_template1 <- raster(mr, xmn = -71, xmx = -55, ymn = 43, ymx = 55, 
    crs = projection(rar)) # create template
plot(rar_template1)


## ----rar_warp2, cache = TRUE------------------------------------------
rar_w1 <- resample(rar, rar_template1)
plot(rar_w1)


## ----rar_warp3, cache = TRUE------------------------------------------
rar_template2 <- raster(xmn = -71, xmx = -55, ymn = 43, ymx = 55, 
  crs = projection(rar), resolution = .25) # create template
plot(resample(rar, rar_template1))


## ----ras_warp, cache = TRUE-------------------------------------------
ras_template1 <- st_as_stars(st_bbox(ras), nx = 21, ny = 21, 
  values = runif(21 * 21)) # create template
ras_w1 <- st_warp(ras, ras_template1)
plot(ras_w1)


## ----ras_warp2, cache = TRUE------------------------------------------
ras_w2 <- st_warp(ras, cellsize = 0.25, crs = st_crs(ras)) 
plot(st_warp(ras, ras_w2))


## ----rasterr1, eval = FALSE-------------------------------------------
## # not run
## stl_r <- rasterize(stl, rar)
## plot(stl_r)


## ----rasters1, cache = TRUE-------------------------------------------
stl_s <- st_rasterize(stl, dy = .1, dx = .1)
plot(stl_s)


## ----rar_stc, cache = TRUE--------------------------------------------
rar_stc <- stack(list(bath_v1 = rar, bath_v2 = rar * 2)) # could be a file path
class(rar_stc) 
nlayers(rar_stc)


## ----rar_stc1b, cache = TRUE------------------------------------------
class(rar_stc[[1]]) 
class(rar_stc[[2]]) 


## ----rar_stc2, cache = TRUE-------------------------------------------
plot(rar_stc)


## ----rar_stc3, cache = TRUE-------------------------------------------
plot(crop(rar_stc, extent(-65, -60, 45, 50)))


## ----rar_stc4, cache = TRUE-------------------------------------------
rar_sum <- calc(rar_stc, fun = sum)
rar_sum


## ----rar_stc5, cache = TRUE-------------------------------------------
plot(rar_sum) 


## ----ras_stc0---------------------------------------------------------
ras


## ----ras_stc1---------------------------------------------------------
ras_stc1 <- c(ras, ras * 2)
names(ras_stc1) <- c("bath_v1", "bath_v2")
ras_stc1


## ----ras_stc1b--------------------------------------------------------
plot(ras_stc1)


## ----ras_stc2---------------------------------------------------------
ras_stc1[1] # or ras_stc1["bath_v1"]


## ----ras_stc2b--------------------------------------------------------
ras_stc1[2] # or ras_stc1["bath_v2"]


## ----ras_stc2c--------------------------------------------------------
library(dplyr)
ras_stc1 %>% select("bath_v1")


## ----ras_stc2d--------------------------------------------------------
ras_stc1 %>% mutate(bath_v3 = bath_v2 * 2)


## ----ras_stc3, cache = TRUE-------------------------------------------
ras_stc2 <- c(ras, ras*2, along = "z")
ras_stc2


## ----ras_stc4, cache = TRUE-------------------------------------------
plot(ras_stc2)


## ----ras_stc5, cache = TRUE-------------------------------------------
ext <- st_bbox(c(xmin = -65, xmax = -60, ymin = 45, ymax = 50), crs = 4326)
plot(st_crop(ras_stc2, ext))


## ----ras_stc6, cache = TRUE-------------------------------------------
(ras_ap <- st_apply(ras_stc2, c(1, 2), sum))


## ----ras_stc6b, cache = TRUE------------------------------------------
plot(ras_ap)


## ----plot0, cache = TRUE----------------------------------------------
plot(ras)


## ---- pre_plot--------------------------------------------------------
# breaks
bks <- c(seq(-5000, 0, 1000), 250, 500, 750, 1000)
# cols 
cls <- c("#c7cbce", "#687677", "#222d3d", "#25364a", "#172434", 
  "#ad6a11", "#e6a331", "#e4be29", "#f2ea8b")


## ---- plot1-----------------------------------------------------------
par(las = 1, bg = "#f79c74", cex.axis = .7, mar = c(2, 2, 2, 4))
plot(ras,  breaks = bks, col = cls, main = "St-Lawrence map", axes = TRUE)


## ---- plot1b----------------------------------------------------------
par(las = 1, bg = "#f79c74", cex.axis = .7, mar = c(2, 2, 2, 4), oma = c(0, 2, 0 ,2))
plot(ras,  breaks = bks, col = cls, main = "St-Lawrence map", axes = TRUE)


## ----tmap0------------------------------------------------------------
library(tmap)
map0 <- tm_shape(stl) + tm_borders(col = "red")
map0


## ----tmap1------------------------------------------------------------
map1 <- map0 + tm_compass(type = "8star", position = c("left", "top")) 
map1 


## ----tmap2, cache = TRUE----------------------------------------------
names(ras)
map2 <- tm_shape(ras) + tm_raster("bathy.tif")


## ----tmap3, cache = TRUE----------------------------------------------
map3 <- map2 + map1  # the order matters
map3


## ----tmap4, cache = TRUE----------------------------------------------
map4 <- map2 + map1 + tm_style("bw")
map4


## ----tmap5, cache = TRUE----------------------------------------------
map5 <- tm_shape(ras) + tm_raster("bathy.tif", breaks = c(seq(-5000,
    0, 1000), 250, 500, 750, 1000), palette = "viridis") 
map5 

