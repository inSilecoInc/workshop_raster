## ----read_raster--------------------------------------------------------------
library(raster)
rar <- raster("data/bathy.tif")  
class(rar)


## ----plot_raster--------------------------------------------------------------
plot(rar)


## ----read_stars---------------------------------------------------------------
library(stars)
ras <- read_stars("data/bathy.tif")   # argument `driver` to specify the driver
class(ras)


## ----plot_stars, cache = TRUE-------------------------------------------------
plot(ras)


## ----raster_drivers-----------------------------------------------------------
ra_dr <- rgdal::gdalDrivers()
head(ra_dr)


## ----raster_drivers2----------------------------------------------------------
ra_dr[which(ra_dr$name == "GTiff"), ]


## ----stars_drivers------------------------------------------------------------
st_dr <- sf::st_drivers(what = "raster")
head(st_dr)


## ----stars_drivers2-----------------------------------------------------------
st_dr[which(st_dr$name == "GTiff"), ]


## ----create_dir---------------------------------------------------------------
dir.create("output", showWarnings = FALSE)


## ----write_raster0------------------------------------------------------------
ra_dr$name[ra_dr$create]


## ----write_raster1, cache = TRUE----------------------------------------------
writeRaster(rar, filename = "output/rar.gpkg", format = "GPKG", overwrite = TRUE)


## ----write_stars0-------------------------------------------------------------
sort(st_dr$name[st_dr$write])


## ----write_stars1, cache = TRUE-----------------------------------------------
write_stars(ras, dsn = "output/ras.gpkg", driver = "GPKG")


## ----sol1aa, cache = TRUE, include = TRUE-------------------------------------
# raster 
rr <- raster("data/bathy.tif")
writeRaster(rr, filename = "output/rar.grd", overwrite = TRUE)


## ----sol1ab, cache = TRUE, include = TRUE-------------------------------------
# stars 
rs <- read_stars("data/bathy.tif")
write_stars(rs, dsn = "output/ras.grd", driver = "rraster")


## ----sol1ba, cache = TRUE, include = TRUE-------------------------------------
# raster 
conv_raster_r <- function(input, output, ...) {
  rr <- raster(input)
  writeRaster(rr, filename = output, ...)
}
conv_raster_r("data/bathy.tif", "output/rar2.grd", overwrite = TRUE)


## ----sol1bb, cache = TRUE, include = TRUE-------------------------------------
# stars
conv_raster_s <- function(from, to, ...) {
  rs <- read_stars("data/bathy.tif")
  write_stars(rs, dsn = "output/ras.grd",...)
}
conv_raster_s("data/bathy.tif", "output/ras2.grd", driver = "rraster")


## ----sol1c, cache = TRUE, include = TRUE--------------------------------------
conv_raster_stars <- function(x) {
  if (class(x) == "RasterLayer") {
    tmp <- tempfile(fileext = ".tif")
    writeRaster(x, tmp)
    read_stars(tmp)
  } else if (class(x) == "stars") {
    tmp <- tempfile(fileext = ".tif")
    write_stars(x, tmp)
    raster(tmp)
  } else {
    warning("x should be of class `stars` or `RasterLayer`")
    NULL
  }
}
ras_n <- conv_raster_stars(rar)
rar_n <- conv_raster_stars(ras)
# conv_raster_stars(1)


## ----raster_obj---------------------------------------------------------------
class(rar)
rar


## ----raster_obj_val-----------------------------------------------------------
values(rar)


## ----raster_obj_fun-----------------------------------------------------------
dim(rar)
ncell(rar)
projection(rar)
bbox(rar)
extent(rar)


## ----stars_obj----------------------------------------------------------------
class(ras)
ras


## ----stars_val----------------------------------------------------------------
class(ras[[1]])
ras[[1]]


## ----stars_obj_fun1-----------------------------------------------------------
dim(ras)
ncell(ras)
st_bbox(ras)


## ----stars_obj_fun2-----------------------------------------------------------
ras_crs <- st_crs(ras)
class(ras_crs)
ras_crs


## ----stars_obj_fun3-----------------------------------------------------------
ras_crs$input
ras_crs$proj4string
ras_crs$epsg


## ----convert1-----------------------------------------------------------------
rar_c <- st_as_stars(rar) 
class(rar_c)


## ----convert2-----------------------------------------------------------------
ras_c <- as(ras, "Raster") 
class(ras_c)


## ----bc_rar1------------------------------------------------------------------
class(rar[,])
rar[1, 1]
rar[1:10, 11:20]
mean(rar[,]) # or mean(as.matrix(rar))
quantile(rar[,])


## ----bc_rar1_2----------------------------------------------------------------
rar2 <- rar # create a copy 
rar2[rar > 0] <- NA  # filter
plot(rar2)


## ----bc_ras1------------------------------------------------------------------
# extract value
class(ras[[1]])
ras[[1]][1, 1]
ras[[1]][1:10, 11:20]


## ----bc_ras1_2----------------------------------------------------------------
mean(ras[[1]]) 
quantile(ras[[1]])


## ----bc_rar1_3----------------------------------------------------------------
ras2 <- ras # create a copy
ras2[[1]][ras[[1]] > units::as_units(0, "m")] <- NA # filter
plot(ras2)


## ----sol2aa, include = TRUE---------------------------------------------------
# raster 
rr_d <- rr_e <- rar 
rr_d[rar > 0] <- NA
rr_e[rar < 0] <- NA
mean(values(rr_d), na.rm = TRUE)
mean(values(rr_e), na.rm = TRUE)


## ----sol2ab, include = TRUE---------------------------------------------------
# stars
rs_e <- rs_d <- ras
rs_d[[1]][ras[[1]] > units::as_units(0, "m")] <- NA # filter
rs_e[[1]][ras[[1]] < units::as_units(0, "m")] <- NA 
mean(rs_d[[1]], na.rm = TRUE)
mean(rs_e[[1]], na.rm = TRUE)


## ----sol2ba1, cache = TRUE, include = TRUE------------------------------------
# raster
mean_long_lat_r <- function(x, gt_na = NULL, lt_na = NULL) {
  m <- as.matrix(x)
  if (!is.null(lt_na)) m[m < lt_na] <- NA
  if (!is.null(gt_na)) m[m > gt_na] <- NA 
  list(
    m_lat = apply(m, 1, mean, na.rm = TRUE),
    m_lon = apply(m, 2, mean, na.rm = TRUE)
  )
} 


## ----sol2ba2, cache = TRUE, fig.width = 7, include = TRUE---------------------
res1 <- mean_long_lat_r(rar)
par(mfrow = c(1, 2))
plot(res1$m_lat, main = "Along latitude", ylab = "elevation", xlab = "lat_index")
plot(res1$m_lon, main = "Along longitude", ylab = "elevation", xlab = "lon_index")    


## ----sol2ba3, cache = TRUE, fig.width = 7, include = TRUE---------------------
res2 <- mean_long_lat_r(rar, gt_na = 0)
par(mfrow = c(1, 2))
plot(res2$m_lat, main = "Along latitude", ylab = "elevation", xlab = "lat_index")
plot(res2$m_lon, main = "Along longitude", ylab = "elevation", xlab = "lon_index")


## ----sol2bb1, include = TRUE--------------------------------------------------
# stars
mean_long_lat_s <- function(x, gt_na = NULL, lt_na = NULL) {
  m <- units::drop_units(ras[[1]])
  if (!is.null(lt_na)) m[m < lt_na] <- NA
  if (!is.null(gt_na)) m[m > gt_na] <- NA 
  list(
    m_lon = apply(m, 1, mean, na.rm = TRUE),
    m_lat = apply(m, 2, mean, na.rm = TRUE)
  )
} 


## ----sol2bb2, cache = TRUE, fig.width = 7, include = TRUE---------------------
res1 <- mean_long_lat_s(ras)
par(mfrow = c(1, 2))
plot(res1$m_lon, main = "Along latitude", ylab = "elevation", xlab = "lat_index")
plot(res1$m_lat, main = "Along longitude", ylab = "elevation", xlab = "lon_index")    


## ----sol2bb3, cache = TRUE, fig.width = 7, include = TRUE---------------------
res2 <- mean_long_lat_s(ras, gt_na = 0)
par(mfrow = c(1, 2))
plot(res2$m_lon, main = "Along latitude", ylab = "elevation", xlab = "lat_index")
plot(res2$m_lat, main = "Along longitude", ylab = "elevation", xlab = "lon_index")    


## ----rar_new_crs--------------------------------------------------------------
rar_crs <- CRS("+proj=merc +a=6378137 +b=6378137 +lat_ts=0.0 +lon_0=0.0 +x_0=0.0 +y_0=0 +k=1.0 +units=m +nadgrids=@null +wktext  +no_defs")


## ----projrar, cache = TRUE----------------------------------------------------
rar_t <- projectRaster(rar, crs = rar_crs)
projection(rar_t)
plot(rar)


## ----projrar2, cache = TRUE---------------------------------------------------
ras_t <- st_transform(ras, crs = 3857)
st_crs(ras_t)


## ----crop0--------------------------------------------------------------------
plot(rar)
rect(-65, 45, -60, 50)


## ----crop_rar, cache = TRUE---------------------------------------------------
rar_c <- crop(rar, extent(-65, -60, 45, 50))
plot(rar_c)


## ----crop_ras, cache = TRUE---------------------------------------------------
ext <- st_bbox(c(xmin = -65, xmax = -60, ymin = 45, ymax = 50), crs = 4326)
ras_c <- st_crop(ras, ext)
plot(ras_c)


## ----stl----------------------------------------------------------------------
(stl <- sf::st_read("data/st_laurence.geojson"))


## ----plot_stl, cache = TRUE---------------------------------------------------
plot(sf::st_geometry(stl))


## ----rar_m, cache = TRUE------------------------------------------------------
rar_m <- mask(rar, stl)
plot(rar_m)


## ----rar_m2, cache = TRUE-----------------------------------------------------
rar_mi <- mask(rar, stl, inverse = TRUE)
plot(rar_mi)


## ----ras_m, cache = TRUE------------------------------------------------------
ras_m <- ras[stl]
plot(rar_m)


## ----ras_m1, cache = TRUE-----------------------------------------------------
stl_i <- sf::st_difference(sf::st_as_sfc(st_bbox(ras)), stl)
plot(stl_i, col = 2)


## ----ras_m2, cache = TRUE-----------------------------------------------------
ras_mi <- ras[stl_i]
plot(ras_mi)


## ----sol3aa, cache = TRUE, include = TRUE-------------------------------------
# raster 
rar_m <- mask(rar, stl)
mean(rar_m[[1]], na.rm = TRUE)


## ----sol3ab, cache = TRUE, include = TRUE-------------------------------------
# stars 
ras_m <- ras[stl]
mean(ras_m[[1]], na.rm = TRUE)


## ----sol3b0, cache = TRUE, include = TRUE-------------------------------------
plot(rar)
rect(-71, 43, -65.1, 55, border = 2, lwd = 4) 
rect(-64.9, 43, -55, 55, lwd = 4) 


## ----sol3b1, cache = TRUE, include = TRUE-------------------------------------
# raster 
rar_m <- mask(rar, stl)
rar_river <- crop(rar_m, extent(-71, -65, 43, 55))
rar_gulf <- crop(rar_m, extent(-65, -55, 43, 55))
mean(values(rar_river), na.rm = TRUE)
mean(values(rar_gulf), na.rm = TRUE)


## ----sol3b1b, cache = TRUE, fig.width = 7.5, echo = -1, include = TRUE--------
par(mfrow = c(1, 2), oma = c(0, 2, 0, 3))
plot(rar_river)
plot(rar_gulf)


## ----sol3b1c, cache = TRUE, include = TRUE------------------------------------
mean(values(rar_river), na.rm = TRUE)
mean(values(rar_gulf), na.rm = TRUE)


## ----sol3b2, cache = TRUE, include = TRUE-------------------------------------
# stars 
ras_m <- ras[stl]
ras_river <- st_crop(ras_m, 
    st_bbox(c(xmin = -71, xmax = -65, ymin = 43, ymax = 55), crs = 4326)
)
ras_gulf <- st_crop(ras_m, 
  st_bbox(c(xmin = -65, xmax = -55, ymin = 43, ymax = 55), crs = 4326)
)
mean(ras_river[[1]], na.rm = TRUE)
mean(ras_gulf[[1]], na.rm = TRUE)


## ----rar_warp, cache = TRUE---------------------------------------------------
# create a raster with the given resolution
mr <- matrix(runif(21*21), 21, 21)
rar_template1 <- raster(mr, xmn = -71, xmx = -55, ymn = 43, ymx = 55, 
    crs = projection(rar)) # create template
plot(rar_template1)


## ----rar_warp2, cache = TRUE--------------------------------------------------
rar_w1 <- resample(rar, rar_template1)
plot(rar_w1)


## ----rar_warp3, cache = TRUE--------------------------------------------------
rar_template2 <- raster(xmn = -71, xmx = -55, ymn = 43, ymx = 55, 
  crs = projection(rar), resolution = .25) # create template
plot(resample(rar, rar_template2))


## ----ras_warp, cache = TRUE---------------------------------------------------
ras_template1 <- st_as_stars(st_bbox(ras), nx = 21, ny = 21, 
  values = runif(21 * 21)) # create template
ras_w1 <- st_warp(ras, ras_template1)
plot(ras_w1)


## ----ras_warp2, cache = TRUE--------------------------------------------------
ras_w2 <- st_warp(ras, cellsize = 0.25, crs = st_crs(ras)) 
plot(ras_w2)


## ----rasterr1, eval = FALSE---------------------------------------------------
## # not run
## stl_r <- rasterize(stl, rar)
## plot(stl_r)


## ----rasters1, cache = TRUE---------------------------------------------------
stl_s <- st_rasterize(stl, dy = .1, dx = .1)
plot(stl_s)


## ----rar_stc, cache = TRUE----------------------------------------------------
rar_stc <- stack(list(bath_v1 = rar, bath_v2 = rar * 2)) # could be a file path
class(rar_stc) 
nlayers(rar_stc)


## ----rar_stc1b, cache = TRUE--------------------------------------------------
class(rar_stc[[1]]) 
class(rar_stc[[2]]) 


## ----rar_stc2, cache = TRUE---------------------------------------------------
plot(rar_stc)


## ----rar_stc3, cache = TRUE---------------------------------------------------
plot(crop(rar_stc, extent(-65, -60, 45, 50)))


## ----rar_stc4, cache = TRUE---------------------------------------------------
rar_sum <- calc(rar_stc, fun = sum)
rar_sum


## ----rar_stc5, cache = TRUE---------------------------------------------------
plot(rar_sum) 


## ----ras_stc0-----------------------------------------------------------------
ras


## ----ras_stc1-----------------------------------------------------------------
ras_stc1 <- c(ras, ras * 2)
names(ras_stc1) <- c("bath_v1", "bath_v2")
ras_stc1


## ----ras_stc1b----------------------------------------------------------------
plot(ras_stc1)


## ----ras_stc2-----------------------------------------------------------------
ras_stc1[1] # or ras_stc1["bath_v1"]


## ----ras_stc2b----------------------------------------------------------------
ras_stc1[2] # or ras_stc1["bath_v2"]


## ----ras_stc2c----------------------------------------------------------------
library(dplyr)
ras_stc1 %>% select("bath_v1")


## ----ras_stc2d----------------------------------------------------------------
ras_stc1 %>% mutate(bath_v3 = bath_v2 * 2)


## ----ras_stc3, cache = TRUE---------------------------------------------------
ras_stc2 <- c(ras, ras*2, along = "z")
ras_stc2


## ----ras_stc4, cache = TRUE---------------------------------------------------
plot(ras_stc2)


## ----ras_stc5, cache = TRUE---------------------------------------------------
ext <- st_bbox(c(xmin = -65, xmax = -60, ymin = 45, ymax = 50), crs = 4326)
plot(st_crop(ras_stc2, ext))


## ----ras_stc6, cache = TRUE---------------------------------------------------
(ras_ap <- st_apply(ras_stc2, c(1, 2), sum))


## ----ras_stc6b, cache = TRUE--------------------------------------------------
plot(ras_ap)


## ----plot0, cache = TRUE------------------------------------------------------
plot(ras)


## ---- pre_plot----------------------------------------------------------------
# breaks
bks <- c(seq(-5000, 0, 1000), 250, 500, 750, 1000)
# cols 
cls <- c("#c7cbce", "#687677", "#222d3d", "#25364a", "#172434", 
  "#ad6a11", "#e6a331", "#e4be29", "#f2ea8b")


## ---- plot1-------------------------------------------------------------------
par(las = 1, bg = "#f79c74", cex.axis = .7, mar = c(2, 2, 2, 4))
plot(ras,  breaks = bks, col = cls, main = "St-Lawrence map", axes = TRUE)


## ---- plot1b------------------------------------------------------------------
par(las = 1, bg = "#f79c74", cex.axis = .7, mar = c(2, 2, 2, 4), oma = c(0, 2, 0 ,2))
plot(ras,  breaks = bks, col = cls, main = "St-Lawrence map", axes = TRUE)


## ---- soldd1------------------------------------------------------------------
# load files
gulf_region <- sf::read_sf("data/st_laurence.geojson")
strs <- stars::read_stars("data/bathy.tif")
# create raster mask
template <- strs
template[[1]][] <- 0
gulf_region$val_ras <- 1

rasterized <- stars::st_rasterize(gulf_region["val_ras"], template = template)
range(rasterized[[1]])


## ---- soldd2------------------------------------------------------------------
plot(rasterized)


## ---- soldd3------------------------------------------------------------------
# apply raster mask
strs[rasterized == 0] <- NA
plot(strs)


## ---- eval = FALSE------------------------------------------------------------
## install.packages("mapedit")
## library(mapedit)
## anticosti <- editMap()


## ---- eval = FALSE------------------------------------------------------------
## stl <- sf::st_read("data/st_laurence.geojson")
## png("output/anticosti.png", width = 7, height = 5, units = "in", res = 300)
## plot(st_geometry(stl))
## plot(st_geometry(anticosti), col = 2, add = TRUE)
## dev.off()


## ----tmap0--------------------------------------------------------------------
library(tmap)
map0 <- tm_shape(stl) + tm_borders(col = "red")
map0


## ----tmap1--------------------------------------------------------------------
map1 <- map0 + tm_compass(type = "8star", position = c("left", "top")) 
map1 


## ----tmap2, cache = TRUE------------------------------------------------------
names(ras)
map2 <- tm_shape(ras) + tm_raster("bathy.tif")


## ----tmap3, cache = TRUE------------------------------------------------------
map3 <- map2 + map1  # the order matters
map3


## ----tmap4, cache = TRUE------------------------------------------------------
map4 <- map2 + map1 + tm_style("bw")
map4


## ----tmap5, cache = TRUE------------------------------------------------------
map5 <- tm_shape(ras) + tm_raster("bathy.tif", breaks = c(seq(-5000,
    0, 1000), 250, 500, 750, 1000), palette = "viridis") 
map5 


## ----tmap_exo, echo = FALSE, cache = TRUE, fig.width = 6, fig.height = 5, purl = TRUE----
# manipulation 
ras_v <- ras[stl]
names(ras_v) <- "Elevation"
# colors 
pal <- colorRampPalette(c("#c7cbce", "#687677", "#222d3d", "#25364a", "#172434", "#ad6a11", "#e6a331", "#e4be29", "#f2ea8b"))
# raster 
elv <- tm_shape(ras_v) + 
  tm_raster("Elevation", breaks = seq(-800, 200, 100), palette = pal(11), midpoint = NA) + 
  tm_layout(main.title = "St-Lawrence river & Gulf", main.title.color = "#ad6a11") + 
  tm_xlab("longitude", size = 0.5) + tm_ylab("latitude", size = 0.5) +
  tm_graticules(lwd = .5, col = "#aaaaaa")
  
shp <- tm_shape(stl) + tm_borders(col = "black", lwd = 2)
oth <- tm_compass(type = "8star", position = c("left", "bottom")) +     
  tm_scale_bar(breaks = c(0, 100, 200), text.size = .8) + 
  tm_logo(c("https://www.r-project.org/logo/Rlogo.png"), position = c("right", "top"), height = 3) 
elv + shp + oth


## ----sol4a, include = TRUE, eval = FALSE, purl = TRUE-------------------------
## # manipulation
## ras_v <- ras[stl]
## names(ras_v) <- "Elevation"
## # colors
## pal <- colorRampPalette(c("#c7cbce", "#687677", "#222d3d", "#25364a",
##   "#172434", "#ad6a11", "#e6a331", "#e4be29", "#f2ea8b"))
## # raster
## elv <- tm_shape(ras_v) +
##   tm_raster("Elevation", breaks = seq(-800, 200, 100), palette = pal(11), midpoint = NA) +
##   tm_layout(main.title = "St-Lawrence (River & Gulf)", main.title.color = "#ad6a11") +
##   tm_xlab("longitude", size = 0.5) + tm_ylab("latitude", size = 0.5) +
##   tm_graticules(lwd = .5, col = "#aaaaaa")


## ----sol4b, include = TRUE, eval = FALSE, purl = TRUE-------------------------
## # borders
## shp <- tm_shape(stl) + tm_borders(col = "black", lwd = 2)
## 
## # other elements
## oth <- tm_compass(type = "8star", position = c("left", "bottom")) +
##     tm_scale_bar(breaks = c(0, 100, 200), text.size = .8) +
##     tm_logo("https://www.r-project.org/logo/Rlogo.png", position = c("right", "top"), height = 3)
## 
## elv + shp + oth

