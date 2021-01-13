# :mortar_board: Rasters with R 
![deploy workshop](https://github.com/inSilecoInc/workshop_raster/workflows/deploy%20workshop/badge.svg)

See https://insilecoinc.github.io/workshop_raster/ 



# Installation 

## Packages 

The following R packages are required : 

- [`sf`](https://CRAN.R-project.org/package=sf)
- [`stars`](https://CRAN.R-project.org/package=stars)
- [`raster`](https://CRAN.R-project.org/package=raster)
- [`rgdal`](https://CRAN.R-project.org/package=rgdal)
- [`tmap`](https://CRAN.R-project.org/package=tmap)

[`sf`](https://CRAN.R-project.org/package=sf) might be trickier to install, check out the [README](https://github.com/r-spatial/sf/) for guidelines to install it. Also, on for Ubuntu users, there is a [comprehensive blog post](https://geocompr.github.io/post/2020/installing-r-spatial-ubuntu/) for installing spatial packages. 


## Local environment

Once installed, you can set your local environment, the following commands may prove helpful: 

- `getwd()` to get the working directory
- `setwd()` to set the working directory
- `dir.create()` to create a directory 
- `list.files()` to liste files in a directory

also the following code can be used to download data and unzip them

```R 
download.file("https://github.com/inSilecoInc/workshop_raster/raw/main/data_and_script.zip", destfile = "tmp.zip")
unzip("tmp.zip")
unlink("tmp.zip") # remove temporary file
```
