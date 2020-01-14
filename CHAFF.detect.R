#--------------------------------------------
#
# Filename: CHAFF.detect.R
#
# Usage:    no initial arguements, edit input file path in code
#
# Purpose:  Detection of chaff in operational NEXRAD moments
#
# Output:   A confidence array ('conf.tot.for.CHAFF', value range = 0-1, continuous) and a flag ('flag.tot.for.CHAFF', value =1 yes, 0 no)
#
# Created:  8.5.2019 dserke
#
#--------------------------------------------

#------------------------
# load required libraries
#------------------------
library(FuzzyR)
library(zoo)
library(pracma)
library(RColorBrewer)
require(graphics)
#library(rgdal)                                                                                                      
library(raster)
library(ggplot2)
library(ncdf4)
library(fields)

#-----------------------
# define constants
#-----------------------
thresh.strong.CHAFF     <- 0.60
thresh.weak.CHAFF       <- 0.60
  
num.radials.for.moveavg <- 20

#-----------------------------
# user defines NEXRAD (polar) or X-BAND (cart)
#-----------------------------
user.radar.band         <- readline(prompt="Enter 'NEXRAD' or 'X-BAND':")
print(paste("User input = ", user.radar.band, sep=""))

#-----------------------------
# load user defined radar file ()
#-----------------------------
if (user.radar.band == "NEXRAD") {
  user.radar.name     <- readline(prompt="Enter NEXRAD name [KFDX or KMTX]")
  print(paste("User input = ", user.radar.name, sep=""))
  if (user.radar.name == "KFDX") {
    #  Case from https://www.thedrive.com/the-war-zone/26827/watch-this-chaff-cloud-rapidly-light-up-radar-over-dozens-of-miles-of-new-mexico-airspace
    # tilt1 (lowest)
    d1 <- nc_open("/d1/serke/projects/CHAFF_detect_ATEC/data/NEXRAD/KFDX/ncswp_KFDX_20190305_202732.552_0.5_SUR_.nc", write = FALSE, verbose = FALSE)
  } else if (user.radar.name == "KMTX") {
    d1 <- nc_open("/d1/serke/projects/CHAFF_detect_ATEC/data/NEXRAD/KMTX/ncf_20191202_233146.nc", write = FALSE, verbose = FALSE)
  }
} else if (user.radar == "X-BAND") {
  # Case from 8/8/2019 @ CRTC (use highest tilt)
  d1 <- nc_open("/d1/serke/projects/CHAFF_detect_ATEC/data/CRTC_OP5R/nc/20190808/ncf_20190808_202203.nc", write = FALSE, verbose = FALSE)
} else {
  print("ERROR")
}

# distribute d1 vars
print(paste("The file has", d1$nvars, "variables"))
d1.var.num               <- seq(1, d1$nvars, by=1)
for (i in 1:length(d1.var.num)) {
  d1.nam <- paste("v", d1.var.num[i], sep = "")
  assign(d1.nam, d1$var[[d1.var.num[i]]])
}  # end of for (i in 1:length(d1.var.num))

# load fields
if (d1$nvars > 12 & user.radar.name == "KFDX") {
  vol.start.time                     <- ncvar_get(d1, v1 )
  range.unambig                      <- ncvar_get(d1, v5 )
  latitude                           <- ncvar_get(d1, v6 )
  longitude                          <- ncvar_get(d1, v7 )
  azimuth                            <- ncvar_get(d1, v38)
  elevation                          <- ncvar_get(d1, v39)
  DZ.tilt1                           <- ncvar_get(d1, v42)
  VE.tilt1                           <- ncvar_get(d1, v43)
  ZDR.tilt1                          <- ncvar_get(d1, v44)
  PHI.tilt1                          <- ncvar_get(d1, v45)
  RHO.tilt1                          <- ncvar_get(d1, v46)
  SW.tilt1                           <- ncvar_get(d1, v47)
} else if (d1$nvars == 12) {
  start_time                         <- ncvar_get(d1, v1)
  stop_time                          <- ncvar_get(d1, v2)
  time_bounds                        <- ncvar_get(d1, v3 )
  lat0                               <- ncvar_get(d1, v4 )
  lon0                               <- ncvar_get(d1, v5)
  grid_mapping_0                     <- ncvar_get(d1, v6)
  DZ.tilt1                           <- ncvar_get(d1, v7)
  VE.tilt1                           <- ncvar_get(d1, v8)
  SW.tilt1                           <- ncvar_get(d1, v9)
  ZDR.tilt1                          <- ncvar_get(d1, v10)
  PHI.tilt1                          <- ncvar_get(d1, v11)
  RHO.tilt1                          <- ncvar_get(d1, v12)
} else if (d1$nvars == 1) {
  vol.start.time                     <- ncvar_get(d1, v1 )
} else if (d1$nvars == 15 & user.radar.name == "KMTX") {
  x0                                 <- ncvar_get(d1, v1 )
  y0                                 <- ncvar_get(d1, v2 )
  lat0                               <- ncvar_get(d1, v4 )
  lon0                               <- ncvar_get(d1, v5 )
  z0                                 <- ncvar_get(d1, v5)
  grid_mapping_0                     <- ncvar_get(d1, v6)
  DZ.tilt1                           <- ncvar_get(d1, v7)[,,4]
  VE.tilt1                           <- ncvar_get(d1, v8)[,,4]
  SW.tilt1                           <- ncvar_get(d1, v9)[,,4]
  ZDR.tilt1                          <- ncvar_get(d1, v10)[,,4]
  PHI.tilt1                          <- ncvar_get(d1, v11)[,,4]
  RHO.tilt1                          <- ncvar_get(d1, v12)[,,4]
  DZ_s3                              <- ncvar_get(d1, v13)
  DZ_s5                              <- ncvar_get(d1, v14)
  DZ_s7                              <- ncvar_get(d1, v15)
  
  #DZ.tilt1                           <- DZ.tilt1[,,4]
} else {
  print("Error")
}
#print(paste("V1 has name", v1$name))
nc_close(d1)

#if (user.radar == "NEXRAD") {
#  #tilt2
#  d2 <- nc_open("/d1/serke/projects/CHAFF_detect_ATEC/data/ncswp_KFDX_20190305_203004.912_1.4_SUR_.nc", write = FALSE, verbose = FALSE)
#  print(paste("The file has", d2$nvars, "variables"))
#  d2.var.num               <- seq(1, d2$nvars, by=1)
#  for (i in 1:length(d2.var.num)) {
#    d2.nam <- paste("v", d2.var.num[i], sep = "")
#    assign(d2.nam, d2$var[[d2.var.num[i]]])
#  }
#  if (d2$nvars >= 9) {
#    vol.start.time                     <- ncvar_get(d2, v1 )
#    range.unambig                      <- ncvar_get(d2, v5 )
#    latitude                           <- ncvar_get(d2, v6 )
#    longitude                          <- ncvar_get(d2, v7 )
#    azimuth                            <- ncvar_get(d2, v38)
#    elevation                          <- ncvar_get(d2, v39)
#    DZ.tilt2                           <- ncvar_get(d2, v42)
#    VE.tilt2                           <- ncvar_get(d2, v43)
#    ZDR.tilt2                          <- ncvar_get(d2, v44)
#    PHI.tilt2                          <- ncvar_get(d2, v45)
#    RHO.tilt2                          <- ncvar_get(d2, v46)
#    SW.tilt2                           <- ncvar_get(d2, v47)
#  } else if (d2$nvars == 1) {
#    vol.start.time                     <- ncvar_get(d2, v1 )
#  }
#  #print(paste("V1 has name", v1$name))
#  nc_close(d2)
#} else if (user.radar == "X-BAND") {
#  print("No other files to load for X-BAND")
#} else {
#  print("ERROR")
#}

range                                <- seq(1.0, 298.75, 0.25)
range2                               <- range[1:720]
DZ.tilt1                             <- DZ.tilt1[1:1192,]
DZ.tilt2                             <- DZ.tilt2[1:1192,]

## set up color palatte
#n_palette <- 11
#marker    <- list(color = rev(brewer.pal(n_palette, "Spectral")))

# compute derived fields
#   std(PHI)
std.PHI.tilt1 <- rollapply(PHI.tilt1, 5, sd, fill=0, partial=TRUE)
#std.PHI.tilt2 <- rollapply(PHI.tilt2, 5, sd, fill=0, partial=TRUE)

#   std(DZ)
std.DZ.tilt1  <- rollapply(DZ.tilt1, 5, sd, fill=0, partial=TRUE)
#std.DZ.tilt2  <- rollapply(DZ.tilt2, 5, sd, fill=0, partial=TRUE)

# compute median values within the defined box
# tilt1
ZDR.med.tilt1       <- median(ZDR.tilt1[36:107,     288:453], na.rm=TRUE)
DZ.med.tilt1        <- median(DZ.tilt1[36:107,      288:453], na.rm=TRUE)
std.DZ.med.tilt1    <- median(std.DZ.tilt1[36:107,  288:453], na.rm=TRUE)
RHO.med.tilt1       <- median(RHO.tilt1[36:107,     288:453], na.rm=TRUE)
PHI.med.tilt1       <- median(PHI.tilt1[36:107,     288:453], na.rm=TRUE)
std.PHI.med.tilt1   <- median(std.PHI.tilt1[36:107, 288:453], na.rm=TRUE)
print(ZDR.med.tilt1)
print(DZ.med.tilt1)
print(std.DZ.med.tilt1)
print(RHO.med.tilt1)
print(PHI.med.tilt1)
print(std.PHI.med.tilt1)

# tilt2
ZDR.med.tilt2       <- median(ZDR.tilt2[36:107,     288:453], na.rm=TRUE)
DZ.med.tilt2        <- median(DZ.tilt2[36:107,      288:453], na.rm=TRUE)
std.DZ.med.tilt2    <- median(std.DZ.tilt2[36:107,  288:453], na.rm=TRUE)
RHO.med.tilt2       <- median(RHO.tilt2[36:107,     288:453], na.rm=TRUE)
PHI.med.tilt2       <- median(PHI.tilt2[36:107,     288:453], na.rm=TRUE)
std.PHI.med.tilt2   <- median(std.PHI.tilt2[36:107, 288:453], na.rm=TRUE)
print(ZDR.med.tilt2)
print(DZ.med.tilt2)
print(std.DZ.med.tilt2)
print(RHO.med.tilt2)
print(PHI.med.tilt2)
print(std.PHI.med.tilt2)

#
ind.DZ.tilt1.na                    <- which(is.na(DZ.tilt1))
DZ.tilt1[ind.DZ.tilt1.na]          <- -99
ind.STDDZ.tilt1.na                 <- which(is.na(std.DZ.tilt1))
std.DZ.tilt1[ind.STDDZ.tilt1.na]   <- -99
ind.ZDR.tilt1.na                   <- which(is.na(ZDR.tilt1))
ZDR.tilt1[ind.ZDR.tilt1.na]        <- -99
ind.RHO.tilt1.na                   <- which(is.na(RHO.tilt1))
RHO.tilt1[ind.ZDR.tilt1.na]        <- -99
ind.PHI.tilt1.na                   <- which(is.na(PHI.tilt1))
PHI.tilt1[ind.ZDR.tilt1.na]        <- -99
ind.STDPHI.tilt1.na                <- which(is.na(std.PHI.tilt1))
std.PHI.tilt1[ind.STDPHI.tilt1.na] <- -99
ind.SW.tilt1.na                    <- which(is.na(SW.tilt1))
SW.tilt1[ind.ZDR.tilt1.na]         <- -99

#
fis <- newfis('test')
fis <- addvar(fis, 'input', 'ZDR', c(-10, 10))
fis <- addmf(fis, 'input', 1, 'ZDR', 'gaussmf', c(2.5, 3))
plotmf(fis, "input", 1)

# define membership functions
#   for DZ
CHAFF.DZ.x.map         <- c(-100, 4, 10, 22, 34, 80)
CHAFF.DZ.y.map         <- c(   0, 0,  1,  1,  0,  0)

#   for std(DZ)
CHAFF.STDDZ.x.map      <- c(-100, 0.5,  1, 3.5, 7, 15)
CHAFF.STDDZ.y.map      <- c(   0,   0,  1,   1, 0,  0)

#   for ZDR
CHAFF.ZDR.x.map        <- c(-100, -8, 0, 6, 8, 20)
CHAFF.ZDR.y.map        <- c(   0,  0, 1, 1, 0,  0)

#   for RHO
CHAFF.RHO.x.map        <- c(-100, 0.33, 0.45, 0.70, 0.82, 1.10)
CHAFF.RHO.y.map        <- c(   0,    0,    1,    1,    0,    0)

#   for PHI
CHAFF.PHI.x.map        <- c(-100, 80, 95, 160, 210, 500)
CHAFF.PHI.y.map        <- c(   0,  0,  1,   1,   0,   0)

#   for std(PHI)
CHAFF.STDPHI.x.map     <- c(-100, 20, 80, 120, 155, 400)
CHAFF.STDPHI.y.map     <- c(   0,  0,  1,   1,   0,   0)

#   for SW
CHAFF.SW.x.map         <- c()
CHAFF.SW.y.map         <- c(0, 1, 1, 0)

# compute interest fields
#   NOTE: interp1 converts matrix to array
d.CHAFF.DZ.int         <- interp1(CHAFF.DZ.x.map,     CHAFF.DZ.y.map,     DZ.tilt1,      method="linear")
d.CHAFF.STDDZ.int      <- interp1(CHAFF.STDDZ.x.map,  CHAFF.STDDZ.y.map,  std.DZ.tilt1,  method="linear")
d.CHAFF.ZDR.int        <- interp1(CHAFF.ZDR.x.map,    CHAFF.ZDR.y.map,    ZDR.tilt1,     method="linear")
d.CHAFF.RHO.int        <- interp1(CHAFF.RHO.x.map,    CHAFF.RHO.y.map,    RHO.tilt1,     method="linear")
d.CHAFF.PHI.int        <- interp1(CHAFF.PHI.x.map,    CHAFF.PHI.y.map,    PHI.tilt1,     method="linear")
d.CHAFF.STDPHI.int     <- interp1(CHAFF.STDPHI.x.map, CHAFF.STDPHI.y.map, std.PHI.tilt1, method="linear")
d.CHAFF.SW.int         <- interp1(CHAFF.SW.x.map,     CHAFF.SW.y.map,     SW.tilt1,      method="linear")

# define weights for all available inputs, based on ???
rel.wts.DZ                  <- 5
rel.wts.STDDZ               <- 5
rel.wts.ZDR                 <- 20
rel.wts.RHO                 <- 30
rel.wts.PHI                 <- 30
rel.wts.STDPHI              <- 30
rel.wts.SW                  <- 0

current.wts.tot             <- rel.wts.DZ + rel.wts.STDDZ + rel.wts.ZDR + rel.wts.RHO + rel.wts.PHI + rel.wts.STDPHI + rel.wts.SW 
print(paste(" TOTAL of rel.wts = ", current.wts.tot, sep=""))

# initialize zeroed flag/confidence arrays.  If interest in CHAFF is above an assigned threshold, then the flag field is set to 1.
conf.tot.for.CHAFF.array    <- zeros(dim(DZ.tilt1)[1], dim(DZ.tilt1)[2])
flag.tot.for.CHAFF.array    <- zeros(dim(DZ.tilt1)[1], dim(DZ.tilt1)[2])

# combined confidence value of interests for CHAFF
d.conf.tot.for.CHAFF.array  <- ((d.CHAFF.DZ.int * rel.wts.DZ) + (d.CHAFF.STDDZ.int * rel.wts.STDDZ) + (d.CHAFF.ZDR.int * rel.wts.ZDR) + (d.CHAFF.RHO.int * rel.wts.RHO) + (d.CHAFF.PHI.int * rel.wts.PHI) + (d.CHAFF.STDPHI.int * rel.wts.STDPHI)) / current.wts.tot
  
# convert array to m x n matrix
d.conf.tot.for.CHAFF.matrix <- matrix(d.conf.tot.for.CHAFF.array, dim(DZ.tilt1)[1], dim(DZ.tilt1)[2])

lat.array                   <-seq(from=min(lat0), to=max(lat0), by=(max(lat0)-min(lat0))/399)
lon.array                   <-seq(from=min(lon0), to=max(lon0), by=(max(lon0)-min(lon0))/399)

#-----------------------------
# plotting ...
#-----------------------------
#   tilt1
par(mfrow=c(2, 3))
par(mar=c(4, 4, 4, 4))
image.plot(lon.array, lat.array, DZ.tilt1, xlim=c(-114, -111), ylim=c(40.0, 42.5), zlim=c(-20, 80), useRaster=TRUE, main="DZ.tilt1", xlab="Lon [deg]", ylab="Lat [deg]")
grid()
text(0.50, 0.50, "*",    cex=1.5)
text(0.55, 0.55, "KMTX", cex=1.5)
image.plot(lon.array, lat.array, std.DZ.tilt1,  xlim=c(-114, -111), ylim=c(40.0, 42.5), zlim=c(0, 10), useRaster=TRUE, main="std.DZ.tilt1", xlab="Lon [deg]", ylab="Lat [deg]")
grid()
text(0.50, 0.50, "*",    cex=1.5)
text(0.55, 0.55, "KMTX", cex=1.5)
image.plot(lon.array, lat.array, ZDR.tilt1, xlim=c(-114, -111), ylim=c(40.0, 42.5), zlim=c(-10, 7),  useRaster=TRUE, main="ZDR.tilt1", xlab="Lon [deg]", ylab="Lat [deg]")
grid()
text(0.50, 0.50, "*",    cex=1.5)
text(0.55, 0.55, "KMTX", cex=1.5)
image.plot(lon.array, lat.array, PHI.tilt1, xlim=c(-114, -111), ylim=c(40.0, 42.5), zlim=c(0, 200), useRaster=TRUE, main="PHI.tilt1", xlab="Lon [deg]", ylab="Lat [deg]")
grid()
text(0.50, 0.50, "*",    cex=1.5)
text(0.55, 0.55, "KMTX", cex=1.5)
image.plot(lon.array, lat.array, std.PHI.tilt1, xlim=c(-114, -111), ylim=c(40.0, 42.5), zlim=c(0, 200), useRaster=TRUE, main="std.PHI.tilt1", xlab="Lon [deg]", ylab="Lat [deg]")
grid()
text(0.50, 0.50, "*",    cex=1.5)
text(0.55, 0.55, "KMTX", cex=1.5)
image.plot(lon.array, lat.array, RHO.tilt1, xlim=c(-114, -111), ylim=c(40.0, 42.5), zlim=c(0, 1),   useRaster=TRUE, main="RHO.tilt1", xlab="Lon [deg]", ylab="Lat [deg]")
grid()
text(0.50, 0.50, "*",    cex=1.5)
text(0.55, 0.55, "KMTX", cex=1.5)

#   confidence in Chaff
par(mfrow=c(1,1))
par(mar=c(4, 4, 4, 4))
image.plot(lon.array, lat.array, d.conf.tot.for.CHAFF.matrix, main="2/12/2019 @ 23:31 UTC, Chaff conf., [algo version 1.0]")
grid()
text(0.50, 0.50, "*",    cex=1.5)
text(0.55, 0.55, "KMTX", cex=1.5)

##   tilt2
#par(mfrow=c(2,2))
#par(mar=c(4, 4, 4, 4))
#image(t(DZ.tilt2),  ylim=c(0, 0.2), zlim=c(-20, 80), useRaster=TRUE, main="DZ.tilt1", xlab="azimuth", ylab="range")
#segments(0.40, 0.03, 0.40, 0.09, col="black")
#segments(0.63, 0.03, 0.63, 0.09, col="black")
#segments(0.40, 0.09, 0.63, 0.09, col="black")
#segments(0.63, 0.03, 0.40, 0.03, col="black")
#image(t(ZDR.tilt2), ylim=c(0, 0.2), zlim=c(-10, 7),  useRaster=TRUE, main="ZDR.tilt1", xlab="azimuth", ylab="range")
#segments(0.40, 0.03, 0.40, 0.09, col="black")
#segments(0.63, 0.03, 0.63, 0.09, col="black")
#segments(0.63, 0.09, 0.40, 0.09, col="black")
#segments(0.63, 0.03, 0.40, 0.03, col="black")
#image(t(PHI.tilt2), ylim=c(0, 0.2), zlim=c(0, 200), useRaster=TRUE, main="PHI.tilt1", xlab="azimuth", ylab="range")
#segments(0.40, 0.03, 0.40, 0.09, col="black")
#segments(0.63, 0.03, 0.63, 0.09, col="black")
#segments(0.63, 0.09, 0.40, 0.09, col="black")
#segments(0.63, 0.03, 0.40, 0.03, col="black")
#image(t(RHO.tilt2), ylim=c(0, 0.2), zlim=c(0, 1),   useRaster=TRUE, main="RHO.tilt1", xlab="azimuth", ylab="range")
#segments(0.40, 0.03, 0.40, 0.09, col="black")
#segments(0.63, 0.03, 0.63, 0.09, col="black")
#segments(0.63, 0.09, 0.40, 0.09, col="black")
#segments(0.63, 0.03, 0.40, 0.03, col="black")
