rm()
dev.off()

library(pracma)
library(ncdf4)
library(sp)
library(raster)
library(lattice)
library(rworldmap)
library(jpeg)
library(fields)
library(spatstat)
library(Matrix)
library(matrixStats)

setwd('C:/Users/baar/Documents/var_var/__manuscript/revision_00_code/get_dem')

## settings

## observable
ncfname = '../input_data/ERA5_15sep2020.nc'
ncin <- nc_open(ncfname)
era5 = list(lon0 = ncvar_get(ncin,"longitude"),lat0 = ncvar_get(ncin,"latitude"))
era5$lat = flipdim(matrix(era5$lat0),dim=1)
era5$obs = flipdim(era5$obs0,dim=2)

## map
worldmap <- getMap(resolution = "high")

## DEM
demname = '../../dem/gtopo30_gis_1km.gri'
rawDEM <- stack(demname)

# extract DEM
dem = list(lon=era5$lon,lat=era5$lat)
Nlon = length(dem$lon)
Nlat = length(dem$lat)
coord = meshgrid(dem$lon,dem$lat)
LON = matrix(coord$X,Nlon*Nlat,1)
LAT = matrix(coord$Y,Nlon*Nlat,1)
ALT_BLUR = matrix(extract(rawDEM$alt_blur,cbind(LON,LAT)),Nlon*Nlat,1)
dem$alt_blur = t(matrix(ALT_BLUR,Nlat,Nlon))

# plot DEM
png(file="dem_alt_blur.png")
filled.contour(dem$lon,dem$lat,dem$alt_blur,plot.axes={plot(worldmap,fill = F, border = "blue",xlim = c(min(era5$lon),max(era5$lon)), ylim = c(min(era5$lat),max(era5$lat)),bg = NA,asp = 1,add=T)})
title('Altitude blur [m]')
dev.off()

# write DEM to disk
saveRDS(dem, file = "../input_data/subsampled_dem.rds")
