# Copyright KNMI 2021
# MIT License file attached
# jouke.de.baar@knmi.nl

rm()
dev.off()

library(pracma)
library(ncdf4)
library(rworldmap)
library(Matrix)
library(geosphere)

setwd('C:/Users/baar/Documents/var_var/__manuscript/revision_00_code')
source('krigFunctions.R')

##############
## settings ##
##############
set.seed(1) # set random seed
Nstation = 128 # 128 # number of stations for sampling
relnoisey = 0.05 # relative noise level for sampling
nxval1 = 256 # 256 # number of brute force iterations for penalty
nxval2 = 256 # 256 # number of brute force iterations for theta

##############
## get data ##
##############

## map
worldmap <- getMap(resolution = "high")

## get observable from ERA5 file
ncfname = 'input_data/ERA5_15sep2020.nc'
ncin <- nc_open(ncfname)
era5 = list(lon0 = ncvar_get(ncin,"longitude"),lat0 = ncvar_get(ncin,"latitude"))
temp = ncvar_get(ncin,"t2m")
era5$obs = apply(temp, c(1,2),mean) # this is where observable is defined
era5$lat = flipdim(matrix(era5$lat),dim=1)
era5$obs = flipdim(era5$obs,dim=2)

## plot observable
png(file="output_figures/era5_obs.png")
filled.contour(era5$lon,era5$lat,era5$obs,plot.axes={plot(worldmap,fill = F, border = "blue",xlim = c(min(era5$lon),max(era5$lon)), ylim = c(min(era5$lat),max(era5$lat)),bg = NA,asp = 1,add=T)})
title('Observed variable (ERA5)')
dev.off()

## DEM
demname = 'input_data/subsampled_dem.rds'
dem = readRDS(demname)

Nlon = length(dem$lon)
Nlat = length(dem$lat)
coord = meshgrid(dem$lon,dem$lat)
dem$LON = matrix(coord$X,Nlon*Nlat,1)
dem$LAT = matrix(coord$Y,Nlon*Nlat,1)
dem$ALT_BLUR = matrix(dem$alt_blur,Nlon*Nlat,1)

png(file="output_figures/dem_alt_blur.png")
filled.contour(dem$lon,dem$lat,dem$alt_blur,plot.axes={plot(worldmap,fill = F, border = "blue",xlim = c(min(era5$lon),max(era5$lon)), ylim = c(min(era5$lat),max(era5$lat)),bg = NA,asp = 1,add=T)})
title('Altitude blur [m]')
dev.off()

##############
## sampling ##
##############

## synthetic stations
station = list(lon=min(era5$lon)+(max(era5$lon)-min(era5$lon))*matrix(rand(Nstation),Nstation,1),lat=min(era5$lat)+(max(era5$lat)-min(era5$lat))*matrix(rand(Nstation),Nstation,1))
station$alt_blur = matrix(interp2(era5$lat,era5$lon,dem$alt_blur,station$lat,station$lon),Nstation,1)
station$obs = matrix(interp2(era5$lat,era5$lon,era5$obs,station$lat,station$lon),Nstation,1)
station$obs = station$obs + relnoisey*std(station$obs)*randn(size(station$obs))

png(file="output_figures/stations.png")
contour(dem$lon,dem$lat,dem$alt_blur,col='grey',xlab='Longitude [deg]',ylab='Latitude [deg]')
plot(worldmap,fill = F, border = "blue",xlim = c(min(era5$lon),max(era5$lon)), ylim = c(min(era5$lat),max(era5$lat)),bg = NA,asp = 1,add=T)
points(station$lon,station$lat,col='red')
dev.off()

x = cbind(station$lon,station$lat)
xi = cbind(dem$LON,dem$LAT)
y = station$obs

################################
## precompute lag for kriging ##
################################
## Algorithm 1 Step 1 and Algorithm 2 Step 1 & Step 2
H = corrLag(x,x,dem,noisey,sigma,c(1,1))

###############################################
## brute force single-objective optimization ##
###############################################
## Algorithm 1 Step 2

## pre-allocate for brute force optimization
pnltyi = linspace(0,30,nxval1)
thetai = exp(linspace(log(1e3),log(200e3),nxval2))
rmseyiStandard = matrix(NA,nxval2,1)
rmseyi = matrix(NA,nxval1,nxval2)
rmsecdfi = matrix(NA,nxval1,nxval2)

print('Standard brute force optimization ...')
for(k in 1:nxval2){
  print(toString(k/nxval2))
  ## Algorithm 1 Step 3
  hyper = c(0.,thetai[k])
  ## Algorithm 1 Step 4 & Step 5
  xval = krigLOOCVfast(x,y,dem,relnoisey,hyper,H$Hhor,H$Hver)
  ## Algorithm 1 Step 6
  rmseyiStandard[k] = xval$rmseyi
}

png(file="output_figures/rmseyiStandard.png")
loglog(thetai,rmseyiStandard,type='b',col='blue',xlab='Correlation length [m]',ylab='x-Validation error [K]')
title('Standard x-validation')
dev.off()

## Algorithm 1 Step 7
idmin = which(rmseyiStandard == min(rmseyiStandard), arr.ind = TRUE)
hyper = c(0.,thetai[idmin])
print(paste('> hyper = ',toString(hyper)))

## compute optimal interpolation
post = krigPost(x,y,dem,relnoisey,hyper,xi)
yiStandard = post$yi
uiStandard = post$ui

png(file="output_figures/yiStandard.png")
filled.contour(dem$lon,dem$lat,t(matrix(yiStandard,Nlat,Nlon)),plot.axes={plot(worldmap,fill = F, border = "blue",xlim = c(min(era5$lon),max(era5$lon)), ylim = c(min(era5$lat),max(era5$lat)),zlim=c(280,320),bg = NA,asp = 1,add=T);points(station$lon,station$lat,col='blue')})
title('Prediction yi')
dev.off()

png(file="output_figures/uiStandard.png")
filled.contour(dem$lon,dem$lat,t(matrix(log10(uiStandard),Nlat,Nlon)),plot.axes={plot(worldmap,fill = F, border = "blue",xlim = c(min(era5$lon),max(era5$lon)), ylim = c(min(era5$lat),max(era5$lat)),bg = NA,asp = 1,add=T);points(station$lon,station$lat,col='blue')})
title('Prediction log10 ui')
dev.off()

##############################################
## brute force multi-objective optimization ##
##############################################
## Algorithm 2 Step 3

print('Proposed brute force optimization ...')
for(j in 1:nxval1){
  print(toString(j/nxval1))
  for(k in 1:nxval2){
    ## Algorithm 2 Step 4
    hyper = c(pnltyi[j],thetai[k])
    ## Algorithm 2 Step 5-8
    xval = krigLOOCVfast(x,y,dem,relnoisey,hyper,H$Hhor,H$Hver)
    ## Algorithm 2 Step 9
    rmseyi[j,k] = xval$rmseyi
    ## Algorithm 2 Step 10
    rmsecdfi[j,k] = xval$rmsecdf
  }
}

## Algorithm 2 Step 11
rmsecdfiCapped = 0*rmseyi + Inf
id = is.finite(rmsecdfi)
rmsecdfiCapped[id] = rmsecdfi[id]
rmse = sqrt( (rmsecdfiCapped/min(rmsecdfiCapped[id]))^2 + (rmseyi/min(rmseyi))^2 )

## Algorithm 2 Step 12
idmin = which(rmse == min(rmse), arr.ind = TRUE)
hyper = c(pnltyi[idmin[1]],thetai[idmin[2]])
print(paste('> hyper = ',toString(hyper)))

## Plot results
png(file="output_figures/rmseyi.png")
filled.contour(pnltyi,log10(thetai),log10(rmseyi),xlab='Penalty [-]',ylab='log10 Correlation length [m]')
title('log10 RMSE yi')
dev.off()

png(file="output_figures/rmseyi_limited.png")
rmseyiLim = rmseyi
limit = 10*min(rmseyi)
rmseyiLim[rmseyiLim>limit] = limit
filled.contour(pnltyi,log10(thetai),log10(rmseyiLim),xlab='Penalty [-]',ylab='log10 Correlation length [m]')
title('log10 RMSE yi (limited)')
dev.off()

png(file="output_figures/rmsecdfi.png")
filled.contour(pnltyi,log10(thetai),log10(rmsecdfi),xlab='Penalty [-]',ylab='log10 Correlation length [m]')
title('log10 RMSE CDF')
dev.off()

png(file="output_figures/rmsecdfi_limited.png")
rmsecdfiLim = rmsecdfi
limit = 10*min(rmsecdfi)
rmsecdfiLim[rmsecdfiLim>limit] = limit
filled.contour(pnltyi,log10(thetai),log10(rmsecdfiLim),xlab='Penalty [-]',ylab='log10 Correlation length [m]')
title('log10 RMSE CDF')
dev.off()

png(file="output_figures/pareto1.png")
id = (rmsecdfi<(3*min(rmsecdfi))) & (rmseyi<(3*min(rmseyi))) 
plot(rmseyi[id],rmsecdfiCapped[id],xlab='RMSE in yi [K]',ylab='RMSE in CDF [K]',col='blue')
title('Pareto plot')
dev.off()

png(file="output_figures/pareto2.png")
filled.contour(pnltyi,log10(thetai),log10(rmse),xlab='Penalty [-]',ylab='log10 Correlation length [m]')
title('log10 RMSE Pareto optimal')
dev.off()

png(file="output_figures/pareto2_limited.png")
rmseLim = rmse
limit = 10*min(rmse)
rmseLim[rmseLim>limit] = limit
filled.contour(pnltyi,log10(thetai),log10(rmseLim),xlab='Penalty [-]',ylab='log10 Correlation length [m]')
title('log10 RMSE Pareto optimal (limited)')
dev.off()


## compute optimal interpolation
post = krigPost(x,y,dem,relnoisey,hyper,xi)
yi = post$yi
ui = post$ui

png(file="output_figures/yi.png")
filled.contour(dem$lon,dem$lat,t(matrix(yi,Nlat,Nlon)),plot.axes={plot(worldmap,fill = F, border = "blue",xlim = c(min(era5$lon),max(era5$lon)), ylim = c(min(era5$lat),max(era5$lat)),zlim=c(280,320),bg = NA,asp = 1,add=T);points(station$lon,station$lat,col='blue')})
title('Prediction yi')
dev.off()
    
png(file="output_figures/ui.png")
filled.contour(dem$lon,dem$lat,t(matrix(log10(ui),Nlat,Nlon)),plot.axes={plot(worldmap,fill = F, border = "blue",xlim = c(min(era5$lon),max(era5$lon)), ylim = c(min(era5$lat),max(era5$lat)),bg = NA,asp = 1,add=T);points(station$lon,station$lat,col='blue')})
title('Prediction log10 ui')
dev.off()
    
