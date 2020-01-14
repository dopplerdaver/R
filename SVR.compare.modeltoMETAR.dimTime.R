#--------------------------------------------
#
# Filename: SVR.compare.modeltoMETAR.dimTime.R
#
# Purposes: Compare metar vis to model vis for case study
#
# Usage:    no initial arguements, code asks six questions of user, need to manually edit input file paths in code
#
# Inputs:   
#
# Process:  
#
# Created:  9.25.2019 dserke
#
# THINGS TO DO:
#           1. would need to output (20 hours) x (4 grids per hour) x (1 height [lowest]) to be able to produce timeseries of point model ~40m vis vs point METAR sfc vis
#           2. run mySQL commands for 10 locations total (2 done, KORD, K3LF)
#
#--------------------------------------------

#------------------------
# load required libraries
#------------------------
#library(FuzzyR)
library(zoo)
library(pracma)
library(RColorBrewer)
require(graphics)
#library(rgdal)   
library(sp)
#library(rspatial) # need to load.package
library(reshape2)
library(raster)
library(ggplot2)
library(ggmap)
library(maps)
#library(mapsdata)
library(ncdf4)
library(NISTunits)
library(lubridate)

require(utils)
require(rgdal)
require(raster)

#-----------------------
# define constants
#-----------------------

VFR.vis.nm              <- 5            # min vis (RVR) [nm] for VFR rules, below which MIFR/IFR
VFR.vis.m               <- NISTmileNauticalTOmeter(VFR.vis.nm)

# case parameters
case.yyyymmdd           <- "20190207"   # yyyymmdd

# model parameters 
NWP.model               <- "WRF"        #
NWP.horiz.res.m         <- 3000         # [m]

# location parameters
lat.KRFD.deg             <-  42.2711    # [deg]
lon.KRFD.deg             <- -89.0940    # [deg]
elev.KRFD.m              <- 226         # [m]

lat.KORD.deg             <-  41.9745219 # [deg]
lon.KORD.deg             <- -87.9065972 # [deg]
elev.KORD.m              <- 207         # [m]
#mod.x.KORD               <- round((5.60/7.42)*1535, digits=0)
#mod.y.KORD               <- round((4.65/9.05)*1023, digits=0)
mod.x.KORD               <- 1031
mod.y.KORD               <- 654

lat.K3LF.deg             <-  39.10      # [deg], Litchfield, IL
lon.K3LF.deg             <- -89.40      # [deg]
elev.K3LF.m              <- 210         # [m]
#mod.x.K3LF               <- round((4.70/7.42)*1535, digits=0)
#mod.y.K3LF               <- round((2.53/9.05)*1023, digits=0)
mod.x.K3LF               <- 1001
mod.y.K3LF               <- 544

lat.KLWA.deg             <-  42.21      # [deg], South Haven, MI
lon.KLWA.deg             <- -86.15      # [deg]
elev.KLWA.m              <- 203         # [m]
mod.x.KLWA               <- 1078 
mod.y.KLWA               <- 669

lat.KDVN.deg             <-  41.37      # [deg], South Haven, MI
lon.KDVN.deg             <- -90.35      # [deg]
elev.KDVN.m              <- 230
mod.x.KDVN               <- 967 
mod.y.KDVN               <- 633         # was 626

lat.KDLH.deg             <-  46.51      # [deg], South Haven, MI
lon.KDLH.deg             <- -92.12      # [deg]
elev.KDLH.m              <- 435
mod.x.KDLH               <- 967 
mod.y.KDLH               <- 633

lat.KAAA.deg             <-  40.10      # [deg], South Haven, MI
lon.KAAA.deg             <- -89.20      # [deg]
elev.KAAA.m              <- 182
mod.x.KAAA               <- 1003 
mod.y.KAAA               <- 582

lat.KMWC.deg             <-  43.07      # [deg], South Haven, MI
lon.KMWC.deg             <- -88.02      # [deg]
elev.KMWC.m              <- 227
mod.x.KMWC               <- 1025 
mod.y.KMWC               <- 695

lat.KDSM.deg             <-  41.32      # [deg], South Haven, MI
lon.KDSM.deg             <- -93.40      # [deg]
elev.KDSM.m              <- 295
mod.x.KDSM               <- 882 
mod.y.KDSM               <- 619

lat.KAZO.deg             <-  42.14      # [deg], South Haven, MI
lon.KAZO.deg             <- -85.33      # [deg]
elev.KAZO.m              <- 272
mod.x.KAZO               <- 1102 
mod.y.KAZO               <- 669

lat.KBDH.deg             <-  45.07+0.2  # [deg], Duluth, MN
lon.KBDH.deg             <- -95.08      # [deg]
elev.KBDH.m              <- 247         # [m]
mod.x.KAZO               <- 832 
mod.y.KAZO               <- 757

hhPERdd                  <- 24
minPERhh                 <- 60
secPERmin                <- 60
days.since.20190207      <- 237

metar.station            <- "KAAA"
mod.x                    <- mod.x.KAAA
mod.y                    <- mod.y.KAAA

#-----------------------------
# convert units
#-----------------------------
mPERkm                  <- 1000

#-----------------------------
# define input data file path and name
#-----------------------------
#   define METAR data from ICICLE domain during case
file.path.metar.visib   <- file.path("/d1/serke/projects/SLANT_VIS_RANGE_DETECT_FAA/data/metar_visib/20190207/")
file.name.metar.visib   <- paste("20190207_", metar.station, "_metar.noraw.txt", sep="")
file.name.metar.visib

#   define WRF-3km vis/lat/lon, lowest-level-only nc file
file.path.wrf.visib     <- file.path("/d1/serke/projects/SLANT_VIS_RANGE_DETECT_FAA/data/wrf_visib/20190207/NUMZ.L1/")
file.name.wrf.T1.visib  <- "20190207_0315.L1.nc"
file.name.wrf.T2.visib  <- "20190207_0400.L1.nc"
file.name.wrf.T3.visib  <- "20190207_0500.L1.nc"
file.name.wrf.T4.visib  <- "20190207_0600.L1.nc"
file.name.wrf.T5.visib  <- "20190207_0700.L1.nc"
file.name.wrf.T6.visib  <- "20190207_0800.L1.nc"
file.name.wrf.T7.visib  <- "20190207_0900.L1.nc"
file.name.wrf.T8.visib  <- "20190207_1000.L1.nc"
file.name.wrf.T9.visib  <- "20190207_1100.L1.nc"
file.name.wrf.T10.visib <- "20190207_1200.L1.nc"
file.name.wrf.T11.visib <- "20190207_1300.L1.nc"
file.name.wrf.T12.visib <- "20190207_1400.L1.nc"
file.name.wrf.T13.visib <- "20190207_1500.L1.nc"
file.name.wrf.T14.visib <- "20190207_1600.L1.nc"
file.name.wrf.T15.visib <- "20190207_1700.L1.nc"
file.name.wrf.T16.visib <- "20190207_1800.L1.nc"
file.name.wrf.T17.visib <- "20190207_1900.L1.nc"
file.name.wrf.T18.visib <- "20190207_2000.L1.nc"
file.name.wrf.T19.visib <- "20190207_2100.L1.nc"
file.name.wrf.T20.visib <- "20190207_2200.L1.nc"

#-----------------------------
# ingest all input data sources
#-----------------------------

# use model level height data from one sample time, a kluge for now but not unreasonable
model.lvl.hgts.atKORD.m <- c(246.76785, 291.20749, 337.75629, 386.51926, 437.59088, 491.11136, 547.19183, 605.91620, 667.39801, 731.76642,  799.16156, 869.70911, 943.59399,  1021.12646,  1102.86487, 1189.14087, 1280.17859,  1375.98535, 1476.66431, 1582.34326,   1693.18872,   1809.44556,   1931.51147,   2059.72827,   2194.46558,   2336.10156,   2485.01196,   2641.55200,   2806.08008,   2978.95239,   3160.50781,   3351.16016,   3551.42725,   3761.85791,   3982.87451,   4214.54248,   4457.07178,   4710.66113,   4975.68408,   5252.60205,   5542.25391,   5845.36621,   6162.43115,   6493.38086,   6823.15625,   7151.84570,   7479.68066,   7806.93213,   8133.75732,   8460.07129,   8785.40820,   9109.31543,   9431.75488,   9754.94531,  10080.73730,  10409.13965,  10739.74805,  11089.48047,  11460.24609,  11853.43555,  12270.32617,  12713.17285,  13185.04980,  13688.41211,  14224.79102,  14796.59766,  15404.41406,  16048.11035,  16731.59766,  17462.07031,  18051.42188)
model.levl.hgts.abvSFC.m<- model.lvl.hgts.atKORD.m - elev.KORD.m

#   ingest WRF-3km vis/lat/lon, lowest-level-only nc file, from case time index #45 (time.model.hhmm[45] = 14:00 UTC)
data.wrf.T1.visib       <- nc_open(paste(file.path.wrf.visib, file.name.wrf.T1.visib, sep=""), write = FALSE, verbose = FALSE)
print(paste("The file has", data.wrf.T1.visib$nvars, "variables"))
data.wrf.T1.visib.var.num<- seq(1, data.wrf.T1.visib$nvars, by=1)
for (i in 1:length(data.wrf.T1.visib.var.num)) {
  data.wrf.T1.visib.nam <- paste("v", data.wrf.T1.visib.var.num[i], sep = "")
  assign(data.wrf.T1.visib.nam, data.wrf.T1.visib$var[[data.wrf.T1.visib.var.num[i]]])
}  # end of for (i in 1:length())
if (data.wrf.T1.visib$nvars == 4) {
  VISIB.wrf.T1.raw             <- ncvar_get( data.wrf.T1.visib, v1 )
  XLAT.wrf.T1.raw              <- ncvar_get( data.wrf.T1.visib, v2 )
  XLON.wrf.T1.raw              <- ncvar_get( data.wrf.T1.visib, v3 )
  XTIME.wrf.T1.raw             <- ncvar_get( data.wrf.T1.visib, v4 )
} else if (data.wrf.T1.visib$nvars == 1) {
  print("Error")
}  # end of if (data.wrf.visib ...)
nc_close(data.wrf.T1.visib)

data.wrf.T2.visib       <- nc_open(paste(file.path.wrf.visib, file.name.wrf.T2.visib, sep=""), write = FALSE, verbose = FALSE)
print(paste("The file has", data.wrf.T2.visib$nvars, "variables"))
data.wrf.T2.visib.var.num<- seq(1, data.wrf.T2.visib$nvars, by=1)
for (i in 1:length(data.wrf.T2.visib.var.num)) {
  data.wrf.T2.visib.nam <- paste("v", data.wrf.T2.visib.var.num[i], sep = "")
  assign(data.wrf.T2.visib.nam, data.wrf.T2.visib$var[[data.wrf.T2.visib.var.num[i]]])
}  # end of for (i in 1:length())
if (data.wrf.T2.visib$nvars == 4) {
  VISIB.wrf.T2.raw             <- ncvar_get( data.wrf.T2.visib, v1 )
  XLAT.wrf.T2.raw              <- ncvar_get( data.wrf.T2.visib, v2 )
  XLON.wrf.T2.raw              <- ncvar_get( data.wrf.T2.visib, v3 )
  XTIME.wrf.T2.raw             <- ncvar_get( data.wrf.T2.visib, v4 )
} else if (data.wrf.T2.visib$nvars == 1) {
  print("Error")
}  # end of if (data.wrf.visib ...)
nc_close(data.wrf.T2.visib)

data.wrf.T3.visib       <- nc_open(paste(file.path.wrf.visib, file.name.wrf.T3.visib, sep=""), write = FALSE, verbose = FALSE)
print(paste("The file has", data.wrf.T3.visib$nvars, "variables"))
data.wrf.T3.visib.var.num<- seq(1, data.wrf.T3.visib$nvars, by=1)
for (i in 1:length(data.wrf.T3.visib.var.num)) {
  data.wrf.T3.visib.nam <- paste("v", data.wrf.T3.visib.var.num[i], sep = "")
  assign(data.wrf.T3.visib.nam, data.wrf.T3.visib$var[[data.wrf.T3.visib.var.num[i]]])
}  # end of for (i in 1:length())
if (data.wrf.T3.visib$nvars == 4) {
  VISIB.wrf.T3.raw             <- ncvar_get( data.wrf.T3.visib, v1 )
  XLAT.wrf.T3.raw              <- ncvar_get( data.wrf.T3.visib, v2 )
  XLON.wrf.T3.raw              <- ncvar_get( data.wrf.T3.visib, v3 )
  XTIME.wrf.T3.raw             <- ncvar_get( data.wrf.T3.visib, v4 )
} else if (data.wrf.T3.visib$nvars == 1) {
  print("Error")
}  # end of if (data.wrf.visib ...)
nc_close(data.wrf.T3.visib)

data.wrf.T4.visib       <- nc_open(paste(file.path.wrf.visib, file.name.wrf.T4.visib, sep=""), write = FALSE, verbose = FALSE)
print(paste("The file has", data.wrf.T4.visib$nvars, "variables"))
data.wrf.T4.visib.var.num<- seq(1, data.wrf.T4.visib$nvars, by=1)
for (i in 1:length(data.wrf.T4.visib.var.num)) {
  data.wrf.T4.visib.nam <- paste("v", data.wrf.T4.visib.var.num[i], sep = "")
  assign(data.wrf.T4.visib.nam, data.wrf.T4.visib$var[[data.wrf.T4.visib.var.num[i]]])
}  # end of for (i in 1:length())
if (data.wrf.T4.visib$nvars == 4) {
  VISIB.wrf.T4.raw             <- ncvar_get( data.wrf.T4.visib, v1 )
  XLAT.wrf.T4.raw              <- ncvar_get( data.wrf.T4.visib, v2 )
  XLON.wrf.T4.raw              <- ncvar_get( data.wrf.T4.visib, v3 )
  XTIME.wrf.T4.raw             <- ncvar_get( data.wrf.T4.visib, v4 )
} else if (data.wrf.T4.visib$nvars == 1) {
  print("Error")
}  # end of if (data.wrf.visib ...)
nc_close(data.wrf.T4.visib)

data.wrf.T5.visib       <- nc_open(paste(file.path.wrf.visib, file.name.wrf.T5.visib, sep=""), write = FALSE, verbose = FALSE)
print(paste("The file has", data.wrf.T5.visib$nvars, "variables"))
data.wrf.T5.visib.var.num<- seq(1, data.wrf.T5.visib$nvars, by=1)
for (i in 1:length(data.wrf.T5.visib.var.num)) {
  data.wrf.T5.visib.nam <- paste("v", data.wrf.T5.visib.var.num[i], sep = "")
  assign(data.wrf.T5.visib.nam, data.wrf.T5.visib$var[[data.wrf.T5.visib.var.num[i]]])
}  # end of for (i in 1:length())
if (data.wrf.T5.visib$nvars == 4) {
  VISIB.wrf.T5.raw             <- ncvar_get( data.wrf.T5.visib, v1 )
  XLAT.wrf.T5.raw              <- ncvar_get( data.wrf.T5.visib, v2 )
  XLON.wrf.T5.raw              <- ncvar_get( data.wrf.T5.visib, v3 )
  XTIME.wrf.T5.raw             <- ncvar_get( data.wrf.T5.visib, v4 )
} else if (data.wrf.T5.visib$nvars == 1) {
  print("Error")
}  # end of if (data.wrf.visib ...)
nc_close(data.wrf.T5.visib)

data.wrf.T6.visib       <- nc_open(paste(file.path.wrf.visib, file.name.wrf.T6.visib, sep=""), write = FALSE, verbose = FALSE)
print(paste("The file has", data.wrf.T6.visib$nvars, "variables"))
data.wrf.T6.visib.var.num<- seq(1, data.wrf.T6.visib$nvars, by=1)
for (i in 1:length(data.wrf.T6.visib.var.num)) {
  data.wrf.T6.visib.nam <- paste("v", data.wrf.T6.visib.var.num[i], sep = "")
  assign(data.wrf.T6.visib.nam, data.wrf.T6.visib$var[[data.wrf.T6.visib.var.num[i]]])
}  # end of for (i in 1:length())
if (data.wrf.T6.visib$nvars == 4) {
  VISIB.wrf.T6.raw             <- ncvar_get( data.wrf.T6.visib, v1 )
  XLAT.wrf.T6.raw              <- ncvar_get( data.wrf.T6.visib, v2 )
  XLON.wrf.T6.raw              <- ncvar_get( data.wrf.T6.visib, v3 )
  XTIME.wrf.T6.raw             <- ncvar_get( data.wrf.T6.visib, v4 )
} else if (data.wrf.T6.visib$nvars == 1) {
  print("Error")
}  # end of if (data.wrf.visib ...)
nc_close(data.wrf.T6.visib)

data.wrf.T7.visib       <- nc_open(paste(file.path.wrf.visib, file.name.wrf.T7.visib, sep=""), write = FALSE, verbose = FALSE)
print(paste("The file has", data.wrf.T7.visib$nvars, "variables"))
data.wrf.T7.visib.var.num<- seq(1, data.wrf.T7.visib$nvars, by=1)
for (i in 1:length(data.wrf.T7.visib.var.num)) {
  data.wrf.T7.visib.nam <- paste("v", data.wrf.T7.visib.var.num[i], sep = "")
  assign(data.wrf.T7.visib.nam, data.wrf.T7.visib$var[[data.wrf.T7.visib.var.num[i]]])
}  # end of for (i in 1:length())
if (data.wrf.T7.visib$nvars == 4) {
  VISIB.wrf.T7.raw             <- ncvar_get( data.wrf.T7.visib, v1 )
  XLAT.wrf.T7.raw              <- ncvar_get( data.wrf.T7.visib, v2 )
  XLON.wrf.T7.raw              <- ncvar_get( data.wrf.T7.visib, v3 )
  XTIME.wrf.T7.raw             <- ncvar_get( data.wrf.T7.visib, v4 )
} else if (data.wrf.T7.visib$nvars == 1) {
  print("Error")
}  # end of if (data.wrf.visib ...)
nc_close(data.wrf.T7.visib)

data.wrf.T8.visib       <- nc_open(paste(file.path.wrf.visib, file.name.wrf.T8.visib, sep=""), write = FALSE, verbose = FALSE)
print(paste("The file has", data.wrf.T8.visib$nvars, "variables"))
data.wrf.T8.visib.var.num<- seq(1, data.wrf.T8.visib$nvars, by=1)
for (i in 1:length(data.wrf.T8.visib.var.num)) {
  data.wrf.T8.visib.nam <- paste("v", data.wrf.T8.visib.var.num[i], sep = "")
  assign(data.wrf.T8.visib.nam, data.wrf.T8.visib$var[[data.wrf.T8.visib.var.num[i]]])
}  # end of for (i in 1:length())
if (data.wrf.T8.visib$nvars == 4) {
  VISIB.wrf.T8.raw             <- ncvar_get( data.wrf.T8.visib, v1 )
  XLAT.wrf.T8.raw              <- ncvar_get( data.wrf.T8.visib, v2 )
  XLON.wrf.T8.raw              <- ncvar_get( data.wrf.T8.visib, v3 )
  XTIME.wrf.T8.raw             <- ncvar_get( data.wrf.T8.visib, v4 )
} else if (data.wrf.T8.visib$nvars == 1) {
  print("Error")
}  # end of if (data.wrf.visib ...)
nc_close(data.wrf.T8.visib)

data.wrf.T9.visib       <- nc_open(paste(file.path.wrf.visib, file.name.wrf.T9.visib, sep=""), write = FALSE, verbose = FALSE)
print(paste("The file has", data.wrf.T9.visib$nvars, "variables"))
data.wrf.T9.visib.var.num<- seq(1, data.wrf.T9.visib$nvars, by=1)
for (i in 1:length(data.wrf.T9.visib.var.num)) {
  data.wrf.T9.visib.nam <- paste("v", data.wrf.T9.visib.var.num[i], sep = "")
  assign(data.wrf.T9.visib.nam, data.wrf.T9.visib$var[[data.wrf.T9.visib.var.num[i]]])
}  # end of for (i in 1:length())
if (data.wrf.T9.visib$nvars == 4) {
  VISIB.wrf.T9.raw             <- ncvar_get( data.wrf.T9.visib, v1 )
  XLAT.wrf.T9.raw              <- ncvar_get( data.wrf.T9.visib, v2 )
  XLON.wrf.T9.raw              <- ncvar_get( data.wrf.T9.visib, v3 )
  XTIME.wrf.T9.raw             <- ncvar_get( data.wrf.T9.visib, v4 )
} else if (data.wrf.T9.visib$nvars == 1) {
  print("Error")
}  # end of if (data.wrf.visib ...)
nc_close(data.wrf.T9.visib)

data.wrf.T10.visib       <- nc_open(paste(file.path.wrf.visib, file.name.wrf.T10.visib, sep=""), write = FALSE, verbose = FALSE)
print(paste("The file has", data.wrf.T10.visib$nvars, "variables"))
data.wrf.T10.visib.var.num<- seq(1, data.wrf.T10.visib$nvars, by=1)
for (i in 1:length(data.wrf.T10.visib.var.num)) {
  data.wrf.T10.visib.nam <- paste("v", data.wrf.T10.visib.var.num[i], sep = "")
  assign(data.wrf.T10.visib.nam, data.wrf.T10.visib$var[[data.wrf.T10.visib.var.num[i]]])
}  # end of for (i in 1:length())
if (data.wrf.T10.visib$nvars == 4) {
  VISIB.wrf.T10.raw             <- ncvar_get( data.wrf.T10.visib, v1 )
  XLAT.wrf.T10.raw              <- ncvar_get( data.wrf.T10.visib, v2 )
  XLON.wrf.T10.raw              <- ncvar_get( data.wrf.T10.visib, v3 )
  XTIME.wrf.T10.raw             <- ncvar_get( data.wrf.T10.visib, v4 )
} else if (data.wrf.T10.visib$nvars == 1) {
  print("Error")
}  # end of if (data.wrf.visib ...)
nc_close(data.wrf.T10.visib)

data.wrf.T11.visib       <- nc_open(paste(file.path.wrf.visib, file.name.wrf.T11.visib, sep=""), write = FALSE, verbose = FALSE)
print(paste("The file has", data.wrf.T11.visib$nvars, "variables"))
data.wrf.T11.visib.var.num<- seq(1, data.wrf.T11.visib$nvars, by=1)
for (i in 1:length(data.wrf.T11.visib.var.num)) {
  data.wrf.T11.visib.nam <- paste("v", data.wrf.T11.visib.var.num[i], sep = "")
  assign(data.wrf.T11.visib.nam, data.wrf.T11.visib$var[[data.wrf.T11.visib.var.num[i]]])
}  # end of for (i in 1:length())
if (data.wrf.T11.visib$nvars == 4) {
  VISIB.wrf.T11.raw             <- ncvar_get( data.wrf.T11.visib, v1 )
  XLAT.wrf.T11.raw              <- ncvar_get( data.wrf.T11.visib, v2 )
  XLON.wrf.T11.raw              <- ncvar_get( data.wrf.T11.visib, v3 )
  XTIME.wrf.T11.raw             <- ncvar_get( data.wrf.T11.visib, v4 )
} else if (data.wrf.T11.visib$nvars == 1) {
  print("Error")
}  # end of if (data.wrf.visib ...)
nc_close(data.wrf.T11.visib)

data.wrf.T12.visib       <- nc_open(paste(file.path.wrf.visib, file.name.wrf.T12.visib, sep=""), write = FALSE, verbose = FALSE)
print(paste("The file has", data.wrf.T12.visib$nvars, "variables"))
data.wrf.T12.visib.var.num<- seq(1, data.wrf.T12.visib$nvars, by=1)
for (i in 1:length(data.wrf.T12.visib.var.num)) {
  data.wrf.T12.visib.nam <- paste("v", data.wrf.T12.visib.var.num[i], sep = "")
  assign(data.wrf.T12.visib.nam, data.wrf.T12.visib$var[[data.wrf.T12.visib.var.num[i]]])
}  # end of for (i in 1:length())
if (data.wrf.T12.visib$nvars == 4) {
  VISIB.wrf.T12.raw             <- ncvar_get( data.wrf.T12.visib, v1 )
  XLAT.wrf.T12.raw              <- ncvar_get( data.wrf.T12.visib, v2 )
  XLON.wrf.T12.raw              <- ncvar_get( data.wrf.T12.visib, v3 )
  XTIME.wrf.T12.raw             <- ncvar_get( data.wrf.T12.visib, v4 )
} else if (data.wrf.T12.visib$nvars == 1) {
  print("Error")
}  # end of if (data.wrf.visib ...)
nc_close(data.wrf.T12.visib)

data.wrf.T13.visib       <- nc_open(paste(file.path.wrf.visib, file.name.wrf.T13.visib, sep=""), write = FALSE, verbose = FALSE)
print(paste("The file has", data.wrf.T13.visib$nvars, "variables"))
data.wrf.T13.visib.var.num<- seq(1, data.wrf.T13.visib$nvars, by=1)
for (i in 1:length(data.wrf.T13.visib.var.num)) {
  data.wrf.T13.visib.nam <- paste("v", data.wrf.T13.visib.var.num[i], sep = "")
  assign(data.wrf.T13.visib.nam, data.wrf.T13.visib$var[[data.wrf.T13.visib.var.num[i]]])
}  # end of for (i in 1:length())
if (data.wrf.T13.visib$nvars == 4) {
  VISIB.wrf.T13.raw             <- ncvar_get( data.wrf.T13.visib, v1 )
  XLAT.wrf.T13.raw              <- ncvar_get( data.wrf.T13.visib, v2 )
  XLON.wrf.T13.raw              <- ncvar_get( data.wrf.T13.visib, v3 )
  XTIME.wrf.T13.raw             <- ncvar_get( data.wrf.T13.visib, v4 )
} else if (data.wrf.T13.visib$nvars == 1) {
  print("Error")
}  # end of if (data.wrf.visib ...)
nc_close(data.wrf.T13.visib)

data.wrf.T14.visib       <- nc_open(paste(file.path.wrf.visib, file.name.wrf.T14.visib, sep=""), write = FALSE, verbose = FALSE)
print(paste("The file has", data.wrf.T14.visib$nvars, "variables"))
data.wrf.T14.visib.var.num<- seq(1, data.wrf.T14.visib$nvars, by=1)
for (i in 1:length(data.wrf.T14.visib.var.num)) {
  data.wrf.T14.visib.nam <- paste("v", data.wrf.T14.visib.var.num[i], sep = "")
  assign(data.wrf.T14.visib.nam, data.wrf.T14.visib$var[[data.wrf.T14.visib.var.num[i]]])
}  # end of for (i in 1:length())
if (data.wrf.T14.visib$nvars == 4) {
  VISIB.wrf.T14.raw             <- ncvar_get( data.wrf.T14.visib, v1 )
  XLAT.wrf.T14.raw              <- ncvar_get( data.wrf.T14.visib, v2 )
  XLON.wrf.T14.raw              <- ncvar_get( data.wrf.T14.visib, v3 )
  XTIME.wrf.T14.raw             <- ncvar_get( data.wrf.T14.visib, v4 )
} else if (data.wrf.T14.visib$nvars == 1) {
  print("Error")
}  # end of if (data.wrf.visib ...)
nc_close(data.wrf.T14.visib)

data.wrf.T15.visib       <- nc_open(paste(file.path.wrf.visib, file.name.wrf.T15.visib, sep=""), write = FALSE, verbose = FALSE)
print(paste("The file has", data.wrf.T15.visib$nvars, "variables"))
data.wrf.T15.visib.var.num<- seq(1, data.wrf.T15.visib$nvars, by=1)
for (i in 1:length(data.wrf.T15.visib.var.num)) {
  data.wrf.T15.visib.nam <- paste("v", data.wrf.T15.visib.var.num[i], sep = "")
  assign(data.wrf.T15.visib.nam, data.wrf.T15.visib$var[[data.wrf.T15.visib.var.num[i]]])
}  # end of for (i in 1:length())
if (data.wrf.T15.visib$nvars == 4) {
  VISIB.wrf.T15.raw             <- ncvar_get( data.wrf.T15.visib, v1 )
  XLAT.wrf.T15.raw              <- ncvar_get( data.wrf.T15.visib, v2 )
  XLON.wrf.T15.raw              <- ncvar_get( data.wrf.T15.visib, v3 )
  XTIME.wrf.T15.raw             <- ncvar_get( data.wrf.T15.visib, v4 )
} else if (data.wrf.T15.visib$nvars == 1) {
  print("Error")
}  # end of if (data.wrf.visib ...)
nc_close(data.wrf.T15.visib)

data.wrf.T16.visib       <- nc_open(paste(file.path.wrf.visib, file.name.wrf.T16.visib, sep=""), write = FALSE, verbose = FALSE)
print(paste("The file has", data.wrf.T16.visib$nvars, "variables"))
data.wrf.T16.visib.var.num<- seq(1, data.wrf.T16.visib$nvars, by=1)
for (i in 1:length(data.wrf.T16.visib.var.num)) {
  data.wrf.T16.visib.nam <- paste("v", data.wrf.T16.visib.var.num[i], sep = "")
  assign(data.wrf.T16.visib.nam, data.wrf.T16.visib$var[[data.wrf.T16.visib.var.num[i]]])
}  # end of for (i in 1:length())
if (data.wrf.T16.visib$nvars == 4) {
  VISIB.wrf.T16.raw             <- ncvar_get( data.wrf.T16.visib, v1 )
  XLAT.wrf.T16.raw              <- ncvar_get( data.wrf.T16.visib, v2 )
  XLON.wrf.T16.raw              <- ncvar_get( data.wrf.T16.visib, v3 )
  XTIME.wrf.T16.raw             <- ncvar_get( data.wrf.T16.visib, v4 )
} else if (data.wrf.T16.visib$nvars == 1) {
  print("Error")
}  # end of if (data.wrf.visib ...)
nc_close(data.wrf.T16.visib)

data.wrf.T17.visib       <- nc_open(paste(file.path.wrf.visib, file.name.wrf.T17.visib, sep=""), write = FALSE, verbose = FALSE)
print(paste("The file has", data.wrf.T17.visib$nvars, "variables"))
data.wrf.T17.visib.var.num<- seq(1, data.wrf.T17.visib$nvars, by=1)
for (i in 1:length(data.wrf.T17.visib.var.num)) {
  data.wrf.T17.visib.nam <- paste("v", data.wrf.T17.visib.var.num[i], sep = "")
  assign(data.wrf.T17.visib.nam, data.wrf.T17.visib$var[[data.wrf.T17.visib.var.num[i]]])
}  # end of for (i in 1:length())
if (data.wrf.T17.visib$nvars == 4) {
  VISIB.wrf.T17.raw             <- ncvar_get( data.wrf.T17.visib, v1 )
  XLAT.wrf.T17.raw              <- ncvar_get( data.wrf.T17.visib, v2 )
  XLON.wrf.T17.raw              <- ncvar_get( data.wrf.T17.visib, v3 )
  XTIME.wrf.T17.raw             <- ncvar_get( data.wrf.T17.visib, v4 )
} else if (data.wrf.T17.visib$nvars == 1) {
  print("Error")
}  # end of if (data.wrf.visib ...)
nc_close(data.wrf.T17.visib)

data.wrf.T18.visib       <- nc_open(paste(file.path.wrf.visib, file.name.wrf.T18.visib, sep=""), write = FALSE, verbose = FALSE)
print(paste("The file has", data.wrf.T18.visib$nvars, "variables"))
data.wrf.T18.visib.var.num<- seq(1, data.wrf.T18.visib$nvars, by=1)
for (i in 1:length(data.wrf.T18.visib.var.num)) {
  data.wrf.T18.visib.nam <- paste("v", data.wrf.T18.visib.var.num[i], sep = "")
  assign(data.wrf.T18.visib.nam, data.wrf.T18.visib$var[[data.wrf.T18.visib.var.num[i]]])
}  # end of for (i in 1:length())
if (data.wrf.T18.visib$nvars == 4) {
  VISIB.wrf.T18.raw             <- ncvar_get( data.wrf.T18.visib, v1 )
  XLAT.wrf.T18.raw              <- ncvar_get( data.wrf.T18.visib, v2 )
  XLON.wrf.T18.raw              <- ncvar_get( data.wrf.T18.visib, v3 )
  XTIME.wrf.T18.raw             <- ncvar_get( data.wrf.T18.visib, v4 )
} else if (data.wrf.T18.visib$nvars == 1) {
  print("Error")
}  # end of if (data.wrf.visib ...)
nc_close(data.wrf.T18.visib)

data.wrf.T19.visib       <- nc_open(paste(file.path.wrf.visib, file.name.wrf.T19.visib, sep=""), write = FALSE, verbose = FALSE)
print(paste("The file has", data.wrf.T19.visib$nvars, "variables"))
data.wrf.T19.visib.var.num<- seq(1, data.wrf.T19.visib$nvars, by=1)
for (i in 1:length(data.wrf.T19.visib.var.num)) {
  data.wrf.T19.visib.nam <- paste("v", data.wrf.T19.visib.var.num[i], sep = "")
  assign(data.wrf.T19.visib.nam, data.wrf.T19.visib$var[[data.wrf.T19.visib.var.num[i]]])
}  # end of for (i in 1:length())
if (data.wrf.T19.visib$nvars == 4) {
  VISIB.wrf.T19.raw             <- ncvar_get( data.wrf.T19.visib, v1 )
  XLAT.wrf.T19.raw              <- ncvar_get( data.wrf.T19.visib, v2 )
  XLON.wrf.T19.raw              <- ncvar_get( data.wrf.T19.visib, v3 )
  XTIME.wrf.T19.raw             <- ncvar_get( data.wrf.T19.visib, v4 )
} else if (data.wrf.T19.visib$nvars == 1) {
  print("Error")
}  # end of if (data.wrf.visib ...)
nc_close(data.wrf.T19.visib)

data.wrf.T20.visib       <- nc_open(paste(file.path.wrf.visib, file.name.wrf.T20.visib, sep=""), write = FALSE, verbose = FALSE)
print(paste("The file has", data.wrf.T20.visib$nvars, "variables"))
data.wrf.T20.visib.var.num<- seq(1, data.wrf.T20.visib$nvars, by=1)
for (i in 1:length(data.wrf.T20.visib.var.num)) {
  data.wrf.T20.visib.nam <- paste("v", data.wrf.T20.visib.var.num[i], sep = "")
  assign(data.wrf.T20.visib.nam, data.wrf.T20.visib$var[[data.wrf.T20.visib.var.num[i]]])
}  # end of for (i in 1:length())
if (data.wrf.T20.visib$nvars == 4) {
  VISIB.wrf.T20.raw             <- ncvar_get( data.wrf.T20.visib, v1 )
  XLAT.wrf.T20.raw              <- ncvar_get( data.wrf.T20.visib, v2 )
  XLON.wrf.T20.raw              <- ncvar_get( data.wrf.T20.visib, v3 )
  XTIME.wrf.T20.raw             <- ncvar_get( data.wrf.T20.visib, v4 )
} else if (data.wrf.T20.visib$nvars == 1) {
  print("Error")
}  # end of if (data.wrf.visib ...)
nc_close(data.wrf.T20.visib)

###########################################################################
## find closest model pixel to point of interest
#ind.sfc.location.lat    <- which.min(abs(XLAT.wrf.T1.raw - lat.K3LF.deg))
#ind.sfc.location.lon    <- which.min(abs(XLON.wrf.T1.raw - lon.K3LF.deg))
#x.num                   <- dim(XLAT.wrf.T1.raw)[1]
#y.num                   <- dim(XLAT.wrf.T1.raw)[2]
#ind.row.num             <- round(ind.sfc.location.lon / x.num)
#ind.col.num             <- round(ind.sfc.location.lat / x.num)

VISIB.wrf.T1.loc            <- VISIB.wrf.T1.raw[mod.x, mod.y]
VISIB.wrf.T2.loc            <- VISIB.wrf.T2.raw[mod.x, mod.y]
VISIB.wrf.T3.loc            <- VISIB.wrf.T3.raw[mod.x, mod.y]
VISIB.wrf.T4.loc            <- VISIB.wrf.T4.raw[mod.x, mod.y]
VISIB.wrf.T5.loc            <- VISIB.wrf.T5.raw[mod.x, mod.y]
VISIB.wrf.T6.loc            <- VISIB.wrf.T6.raw[mod.x, mod.y]
VISIB.wrf.T7.loc            <- VISIB.wrf.T7.raw[mod.x, mod.y]
VISIB.wrf.T8.loc            <- VISIB.wrf.T8.raw[mod.x, mod.y]
VISIB.wrf.T9.loc            <- VISIB.wrf.T9.raw[mod.x, mod.y]
VISIB.wrf.T10.loc           <- VISIB.wrf.T10.raw[mod.x, mod.y]
VISIB.wrf.T11.loc           <- VISIB.wrf.T11.raw[mod.x, mod.y]
VISIB.wrf.T12.loc           <- VISIB.wrf.T12.raw[mod.x, mod.y]
VISIB.wrf.T13.loc           <- VISIB.wrf.T13.raw[mod.x, mod.y]
VISIB.wrf.T14.loc           <- VISIB.wrf.T14.raw[mod.x, mod.y]
VISIB.wrf.T15.loc           <- VISIB.wrf.T15.raw[mod.x, mod.y]
VISIB.wrf.T16.loc           <- VISIB.wrf.T16.raw[mod.x, mod.y]
VISIB.wrf.T17.loc           <- VISIB.wrf.T17.raw[mod.x, mod.y]
VISIB.wrf.T18.loc           <- VISIB.wrf.T18.raw[mod.x, mod.y]
VISIB.wrf.T19.loc           <- VISIB.wrf.T19.raw[mod.x, mod.y]
VISIB.wrf.T20.loc           <- VISIB.wrf.T20.raw[mod.x, mod.y]
visib.m                     <- c(VISIB.wrf.T1.loc, VISIB.wrf.T2.loc, VISIB.wrf.T3.loc, VISIB.wrf.T4.loc, VISIB.wrf.T5.loc, VISIB.wrf.T6.loc, VISIB.wrf.T7.loc, VISIB.wrf.T8.loc, VISIB.wrf.T9.loc, VISIB.wrf.T10.loc, VISIB.wrf.T11.loc, VISIB.wrf.T12.loc, VISIB.wrf.T13.loc, VISIB.wrf.T14.loc, VISIB.wrf.T15.loc, VISIB.wrf.T16.loc, VISIB.wrf.T17.loc, VISIB.wrf.T18.loc, VISIB.wrf.T19.loc, VISIB.wrf.T20.loc)

#
#visib.m                     <- c(  14000,   15000,   14700,   14000,   13000,   12500,   12500,   12300,   13000,   12700,   13000,   13200,   13200,   13800,   13000,   13000,   13100,   12800,   11400,   10900,   10300,    9980,    9900,     500,     490,     300,     500,     400,     340,     400,     420,     400,     400,     500,     700,     400,     290,     250,     250,     190,     180,     170,     200,     200,     180,     150,     160,     170,     190,     210,     280,     400,     600,     800,    1070,    6000,    7800,    7500,    7700,     500,     180,     140,      90,     110,      80,      65,      75,      70,      70,     120,     190,     400,    7000,    8000,   11300,   12500,  13000)
time.model.hhmm             <- c("03:15", "04:00", "05:00", "06:00", "07:00", "08:00", "09:00", "10:00", "11:00", "12:00", "13:00", "14:00", "15:00", "16:00", "17:00", "18:00", "19:00", "20:00", "21:00", "22:00")
#time.model.hhmm             <- c("03:00", "03:15", "03:30", "03:45", "04:00", "04:15", "04:30", "04:45", "05:00", "05:15", "05:30", "05:45", "06:00", "06:15", "06:30", "06:45", "07:00", "07:15", "07:30", "07:45", "08:00", "08:15", "08:30", "08:45", "09:00", "09:15", "09:30", "09:45", "10:00", "10:15", "10:30", "10:45", "11:00", "11:15", "11:30", "11:45", "12:00", "12:15", "12:30", "12:45", "13:00", "13:15", "13:30", "13:45", "14:00", "14:15", "14:30", "14:45", "15:00", "15:15", "15:30", "15:45", "16:00", "16:15", "16:30", "16:45", "17:00", "17:15", "17:30", "17:45", "18:00", "18:15", "18:30", "18:45", "19:00", "19:15", "19:30", "19:45", "20:00", "20:15", "20:30", "20:45", "21:00", "21:15", "21:30", "21:45", "22:00")
time.index                  <- seq(from=1, to=length(visib.m), by=1)
data.model.visib            <- data.frame(time.model.hhmm, time.index, visib.m)
data.model.visib$HH.MM      <- as.POSIXct(strptime(data.model.visib$time.model.hhmm, format="%H:%M"))
data.model.visib$obs        <- replicate(length(data.model.visib$visib.m), "model")
head(data.model.visib)

#   ingest METAR vis obs data file
#     WAITING FOR GREG TO EDIT SCRIPT AND RUN FOR ICICLE DOMAIN
data.metar.visib            <- read.csv(paste(file.path.metar.visib, file.name.metar.visib, sep = ""), sep=" ", header=FALSE)
dim(data.metar.visib)
colnames(data.metar.visib)  <- c("unix.time.s", "obs.date.yyyymmdd", "obs.time.hhmmss", "visib.m", "cld.type.1", "cld.hgt.m.1", "cld.type.2", "cld.hgt.m.2", "cld.type.3", "cld.hgt.m.3")
data.metar.visib$visib.m    <- data.metar.visib$visib.m * 1000
data.metar.visib$HH.MM      <- as.POSIXct(strptime(data.metar.visib$obs.time.hhmmss, format="%H:%M"))
##### UNCOMMENT THIS FOR KBDH
#data.metar.visib$time.index <- c(1, NA, NA, NA, NA, 2, NA, NA, NA, NA, 3, NA, NA, NA, 4, NA, NA, NA, 5, NA, NA, NA, 6, NA, NA, NA, NA, 7, NA, NA, NA, NA, 8, NA, NA, NA, NA, 9, NA, NA, NA, NA, 10, NA, NA, NA, NA, 11, NA, NA, NA, NA, 12, NA, NA, NA, NA, 13, NA, NA, NA, NA, 14, NA, NA, NA, NA, 15, NA, NA, NA, NA, 16, NA, NA, NA, NA, 17, NA, NA, NA, NA, 18, NA, NA, NA, NA, 19, NA, NA, NA, NA, 20)
#num.metar.lines             <- 93
##### UNCOMMENT THIS FOR KAZO
#data.metar.visib$time.index <- c(1, 2, NA, 3, 4, 5, 6, NA, 7, NA, 8, NA, NA, 9, NA, 10, NA, 11, NA, NA, 12, 13, 14, 15, 16, NA, 17, NA, 18, 19, 20)
#num.metar.lines             <- 31
##### UNCOMMENT THIS FOR KDSM
#data.metar.visib$time.index <- c(2, NA, NA, 3, 4, NA, 5, NA, NA, NA, 6, NA, NA, NA, NA, NA, 7, NA, NA, 8, 9, NA, 10, 11, NA, NA, NA, 12, NA, NA, NA, 13, NA, NA, NA, 14, NA, NA, NA, 15, NA, 16, NA, NA, 17, NA, NA, 18, 19, NA, NA, 20)
#num.metar.lines             <- 52
##### UNCOMMENT THIS FOR KMWC
#data.metar.visib$time.index <- c(1, 2, NA, 3, 4, NA, 5, 6, 7, NA, NA, 8, NA, 9, 10, NA, NA, NA, 11, NA, NA, 12, 13, NA, NA, NA, 14, NA, NA, NA, NA, 15, NA, NA, NA, NA, NA, 16, NA, NA, NA, NA, 17, NA, NA, NA, 18, NA, NA, NA, NA, 19, NA, NA, NA, 20, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA)
#num.metar.lines             <- 72
##### UNCOMMENT THIS FOR KAAA
data.metar.visib$time.index <- c(2, NA, 3, NA, NA, 4, NA, NA, 5, NA, NA, 6, NA, NA, 7, NA, NA, 8, NA, NA, 9, NA, NA, 10, NA, NA, 11, NA, NA, 12, NA, NA, 13, NA, NA, 14, NA, NA, 15, NA, NA, 16, NA, NA, 17, NA, NA, 18, NA, NA, 19, NA, NA, 20, NA, NA, NA)
num.metar.lines             <- 57
##### UNCOMMENT THIS FOR KDLH
#data.metar.visib$time.index <- c(NA, NA, 2, NA, 3, NA, NA, NA, 4, NA, NA, NA, NA, 5, NA, 6, NA, NA, 7, NA, NA, 8, NA, NA, 9, NA, NA, NA, 10, NA, 11, NA, NA, 12, NA, NA, NA, NA, 13, NA, NA, NA, 14, NA, NA, NA, 15, NA, 16, NA, NA, NA, 17, 18, 19, NA, 20)
#num.metar.lines             <- 57
##### UNCOMMENT THIS FOR KDVN
#data.metar.visib$time.index <- c(NA, 2, NA, 3, 4, NA, 5, NA, NA, 6, NA, NA, 7, NA, NA, NA, NA, NA, 8, NA, NA, NA, 9, 10, NA, 11, NA, 12, NA, NA, NA, 13, NA, NA, NA, 14, NA, NA, NA, NA, NA, 15, NA, NA, 16, NA, NA, NA, NA, 17, NA, NA, 18, 19, NA, 20)
#num.metar.lines             <- 56
##### UNCOMMENT THIS FOR KLWA
#data.metar.visib$time.index <- c(NA, NA, 2, NA, NA, 3, NA, NA, 4, NA, NA, 5, NA, NA, 6, NA, NA, 7, NA, NA, 8, NA, NA, 9, NA, NA, 10, NA, NA, 11, NA, NA, 12, NA, NA, 13, NA, NA, 14, NA, NA, 15, NA, NA, 16, NA, NA, 17, NA, NA, 18, NA, NA, 19, NA, NA, 20, NA)
#num.metar.lines             <- 58
##### UNCOMMENT THIS FOR K3LF
#data.metar.visib$time.index <- c(1, NA, 2, NA, NA, 3, NA, NA, 4, NA, NA, 5, NA, NA, 6, NA, NA, 7, NA, NA, 8, NA, NA, 9, NA, NA, 10, NA, NA, 11, NA, NA, 12, NA, NA, 13, NA, NA, 14, NA, NA, 15, NA, NA, 16, NA, NA, 17, NA, NA, 18, NA, NA, 19, NA, NA, 20, NA)
#num.metar.lines           <- 58
##### UNCOMMENT THIS FOR KORD DF
#data.metar.visib$time.index <- c(2, NA, NA, 3, NA, NA, 4, NA, 5, NA, NA, 6, 7, NA, NA, NA, 8, NA, NA, 9, 10, NA, NA, NA, 11, NA, 12, NA, NA, 13, 14, 15, NA, NA, 16, 17, NA, 18, NA, NA, NA, NA, NA, 19, NA, NA, 20, NA, NA, NA, NA, NA, NA )
#data.metar.visib$time.index <- c(5, 7, 9, 10, 11, 12, 16, 19, 22, 23, 24, 25, 26, 27, 28, 29, 29.2, 30, 32, 32.5, 36, 39, 40, 40.3, 41, 44, 45, 45.5, 47, 48, 50, 54, 58, 58.3, 59, 62, 65, 66, 68, 68.3, 70, 70.3, 71, 71.2, 73, 74, 75, 79, 83, 87, 91, 95, 99)
#num.metar.lines           <- 53

data.metar.visib$obs        <- replicate(length(data.metar.visib$visib.m), "metar")
head(data.metar.visib)

# make df of rvr values from raw metars
# UNCOMMENT THIS FOR KORD DF
#time.index                  <- c(     20,    20.3,      24,      26,      28,    28.3,    28.6,    29.5,    49.6,    52.4,    56.4,      60,    60.5,    61.5,    64.4,      67)
#HH.MM                       <- c("07:43", "07:51", "08:51", "09:15", "09:44", "09:51", "09:53", "10:07", "15:08", "15:51", "16:51", "17:46", "17:51", "18:07", "18:51", "19:30")
#HH.MM                       <- as.POSIXct(strptime(HH.MM, format="%H:%M"))
#visib.ft                    <- c(   6000,    6000,    5000,    5500,    6000,    6000,    6000,    5500,    6000,    4000,    5000,    6000,    6000,    3000,    6000,    5000)
#visib.m                     <- visib.ft * 0.3048
#data.rvr.visib              <- data.frame(time.index, HH.MM, visib.m)
#data.rvr.visib$obs          <- replicate(length(data.rvr.visib$visib.m), "rvr")
#head(data.rvr.visib)

# make combined df of model, metar-vis and metar-rvr
combine.time.index          <- c(data.model.visib$time.index, data.metar.visib$time.index[1:num.metar.lines])
combine.HH.MM               <- c(data.model.visib$HH.MM,      data.metar.visib$HH.MM[1:num.metar.lines]) - hhPERdd*minPERhh*secPERmin*days.since.20190207 + minPERhh*secPERmin
combine.obs                 <- c(data.model.visib$obs,        data.metar.visib$obs[1:num.metar.lines])
combine.visib               <- c(data.model.visib$visib.m,    data.metar.visib$visib.m[1:num.metar.lines])
data.combine.visib          <- data.frame(combine.time.index, combine.HH.MM, combine.obs, combine.visib )
head(data.combine.visib)
## below is for KORD, where the above data.rvr.visib df is available
#combine.time.index          <- c(data.model.visib$time.index, data.metar.visib$time.index[1:48], data.rvr.visib$time.index)
#combine.HH.MM               <- c(data.model.visib$HH.MM,      data.metar.visib$HH.MM[1:48],      data.rvr.visib$HH.MM) - hhPERdd*minPERhh*secPERmin*days.since.20190207 + minPERhh*secPERmin
#combine.obs                 <- c(data.model.visib$obs,        data.metar.visib$obs[1:48],        data.rvr.visib$obs)
#combine.visib               <- c(data.model.visib$visib.m,    data.metar.visib$visib.m[1:48],    data.rvr.visib$visib.m)
#data.combine.visib          <- data.frame(combine.time.index, combine.HH.MM, combine.obs, combine.visib )

#-----------------------------
# restrict VISIB.raw to less than VFR RVR minimum [m], mainly to restrict the dynamic range of subsequent plots
#-----------------------------
ind.T1.raw.gtVFR                             <- which(VISIB.wrf.T1.raw > VFR.vis.m)
VISIB.wrf.T1.raw[ind.T1.raw.gtVFR]           <- VFR.vis.m
ind.T2.raw.gtVFR                             <- which(VISIB.wrf.T2.raw > VFR.vis.m)
VISIB.wrf.T2.raw[ind.T2.raw.gtVFR]           <- VFR.vis.m
ind.T3.raw.gtVFR                             <- which(VISIB.wrf.T3.raw > VFR.vis.m)
VISIB.wrf.T3.raw[ind.T3.raw.gtVFR]           <- VFR.vis.m
ind.T4.raw.gtVFR                             <- which(VISIB.wrf.T4.raw > VFR.vis.m)
VISIB.wrf.T4.raw[ind.T4.raw.gtVFR]           <- VFR.vis.m
ind.T5.raw.gtVFR                             <- which(VISIB.wrf.T5.raw > VFR.vis.m)
VISIB.wrf.T5.raw[ind.T5.raw.gtVFR]           <- VFR.vis.m
ind.T6.raw.gtVFR                             <- which(VISIB.wrf.T6.raw > VFR.vis.m)
VISIB.wrf.T6.raw[ind.T6.raw.gtVFR]           <- VFR.vis.m
ind.T7.raw.gtVFR                             <- which(VISIB.wrf.T7.raw > VFR.vis.m)
VISIB.wrf.T7.raw[ind.T7.raw.gtVFR]           <- VFR.vis.m
ind.T8.raw.gtVFR                             <- which(VISIB.wrf.T8.raw > VFR.vis.m)
VISIB.wrf.T8.raw[ind.T8.raw.gtVFR]           <- VFR.vis.m
ind.T9.raw.gtVFR                             <- which(VISIB.wrf.T9.raw > VFR.vis.m)
VISIB.wrf.T9.raw[ind.T9.raw.gtVFR]           <- VFR.vis.m
ind.T10.raw.gtVFR                             <- which(VISIB.wrf.T10.raw > VFR.vis.m)
VISIB.wrf.T10.raw[ind.T10.raw.gtVFR]          <- VFR.vis.m
ind.T11.raw.gtVFR                             <- which(VISIB.wrf.T11.raw > VFR.vis.m)
VISIB.wrf.T11.raw[ind.T11.raw.gtVFR]          <- VFR.vis.m
ind.T12.raw.gtVFR                             <- which(VISIB.wrf.T12.raw > VFR.vis.m)
VISIB.wrf.T12.raw[ind.T12.raw.gtVFR]          <- VFR.vis.m
ind.T11.raw.gtVFR                             <- which(VISIB.wrf.T11.raw > VFR.vis.m)
VISIB.wrf.T11.raw[ind.T11.raw.gtVFR]          <- VFR.vis.m
ind.T12.raw.gtVFR                             <- which(VISIB.wrf.T12.raw > VFR.vis.m)
VISIB.wrf.T12.raw[ind.T12.raw.gtVFR]          <- VFR.vis.m
ind.T13.raw.gtVFR                             <- which(VISIB.wrf.T13.raw > VFR.vis.m)
VISIB.wrf.T13.raw[ind.T13.raw.gtVFR]          <- VFR.vis.m
ind.T14.raw.gtVFR                             <- which(VISIB.wrf.T14.raw > VFR.vis.m)
VISIB.wrf.T14.raw[ind.T14.raw.gtVFR]          <- VFR.vis.m
ind.T15.raw.gtVFR                             <- which(VISIB.wrf.T15.raw > VFR.vis.m)
VISIB.wrf.T15.raw[ind.T15.raw.gtVFR]          <- VFR.vis.m
ind.T16.raw.gtVFR                             <- which(VISIB.wrf.T16.raw > VFR.vis.m)
VISIB.wrf.T16.raw[ind.T16.raw.gtVFR]          <- VFR.vis.m
ind.T17.raw.gtVFR                             <- which(VISIB.wrf.T17.raw > VFR.vis.m)
VISIB.wrf.T17.raw[ind.T17.raw.gtVFR]          <- VFR.vis.m
ind.T18.raw.gtVFR                             <- which(VISIB.wrf.T18.raw > VFR.vis.m)
VISIB.wrf.T18.raw[ind.T18.raw.gtVFR]          <- VFR.vis.m
ind.T19.raw.gtVFR                             <- which(VISIB.wrf.T19.raw > VFR.vis.m)
VISIB.wrf.T19.raw[ind.T19.raw.gtVFR]          <- VFR.vis.m
ind.T20.raw.gtVFR                             <- which(VISIB.wrf.T20.raw > VFR.vis.m)
VISIB.wrf.T20.raw[ind.T20.raw.gtVFR]          <- VFR.vis.m

# create model time# df
model.T1.df                                  <- data.frame(as.vector(VISIB.wrf.T1.raw), as.vector(XLON.wrf.T1.raw), as.vector(XLAT.wrf.T1.raw))
colnames(model.T1.df)                        <- c("visib.m", "long", "lat")
#head(model.T1.df)
model.T2.df                                  <- data.frame(as.vector(VISIB.wrf.T2.raw), as.vector(XLON.wrf.T2.raw), as.vector(XLAT.wrf.T2.raw))
colnames(model.T2.df)                        <- c("visib.m", "long", "lat")
#head(model.T2.df)
model.T3.df                                  <- data.frame(as.vector(VISIB.wrf.T3.raw), as.vector(XLON.wrf.T3.raw), as.vector(XLAT.wrf.T3.raw))
colnames(model.T3.df)                        <- c("visib.m", "long", "lat")
#head(model.T3.df)
model.T4.df                                  <- data.frame(as.vector(VISIB.wrf.T4.raw), as.vector(XLON.wrf.T4.raw), as.vector(XLAT.wrf.T4.raw))
colnames(model.T4.df)                        <- c("visib.m", "long", "lat")
#head(model.T4.df)
model.T5.df                                  <- data.frame(as.vector(VISIB.wrf.T5.raw), as.vector(XLON.wrf.T5.raw), as.vector(XLAT.wrf.T5.raw))
colnames(model.T5.df)                        <- c("visib.m", "long", "lat")
#head(model.T5.df)
model.T6.df                                  <- data.frame(as.vector(VISIB.wrf.T6.raw), as.vector(XLON.wrf.T6.raw), as.vector(XLAT.wrf.T6.raw))
colnames(model.T6.df)                        <- c("visib.m", "long", "lat")
#head(model.T6.df)
model.T7.df                                  <- data.frame(as.vector(VISIB.wrf.T7.raw), as.vector(XLON.wrf.T7.raw), as.vector(XLAT.wrf.T7.raw))
colnames(model.T7.df)                        <- c("visib.m", "long", "lat")
#head(model.T7.df)
model.T8.df                                  <- data.frame(as.vector(VISIB.wrf.T8.raw), as.vector(XLON.wrf.T8.raw), as.vector(XLAT.wrf.T8.raw))
colnames(model.T8.df)                        <- c("visib.m", "long", "lat")
#head(model.T8.df)
model.T9.df                                  <- data.frame(as.vector(VISIB.wrf.T9.raw), as.vector(XLON.wrf.T9.raw), as.vector(XLAT.wrf.T9.raw))
colnames(model.T9.df)                        <- c("visib.m", "long", "lat")
#head(model.T9.df)
model.T10.df                                  <- data.frame(as.vector(VISIB.wrf.T10.raw), as.vector(XLON.wrf.T10.raw), as.vector(XLAT.wrf.T10.raw))
colnames(model.T10.df)                        <- c("visib.m", "long", "lat")
#head(model.T10.df)
model.T11.df                                  <- data.frame(as.vector(VISIB.wrf.T11.raw), as.vector(XLON.wrf.T11.raw), as.vector(XLAT.wrf.T11.raw))
colnames(model.T11.df)                        <- c("visib.m", "long", "lat")
#head(model.T11.df)
model.T12.df                                  <- data.frame(as.vector(VISIB.wrf.T12.raw), as.vector(XLON.wrf.T12.raw), as.vector(XLAT.wrf.T12.raw))
colnames(model.T12.df)                        <- c("visib.m", "long", "lat")
#head(model.T12.df)
model.T13.df                                  <- data.frame(as.vector(VISIB.wrf.T13.raw), as.vector(XLON.wrf.T13.raw), as.vector(XLAT.wrf.T13.raw))
colnames(model.T13.df)                        <- c("visib.m", "long", "lat")
#head(model.T13.df)
model.T14.df                                  <- data.frame(as.vector(VISIB.wrf.T14.raw), as.vector(XLON.wrf.T14.raw), as.vector(XLAT.wrf.T14.raw))
colnames(model.T14.df)                        <- c("visib.m", "long", "lat")
#head(model.T14.df)
model.T15.df                                  <- data.frame(as.vector(VISIB.wrf.T15.raw), as.vector(XLON.wrf.T15.raw), as.vector(XLAT.wrf.T15.raw))
colnames(model.T15.df)                        <- c("visib.m", "long", "lat")
#head(model.T15.df)
model.T16.df                                  <- data.frame(as.vector(VISIB.wrf.T16.raw), as.vector(XLON.wrf.T16.raw), as.vector(XLAT.wrf.T16.raw))
colnames(model.T16.df)                        <- c("visib.m", "long", "lat")
#head(model.T16.df)
model.T17.df                                  <- data.frame(as.vector(VISIB.wrf.T17.raw), as.vector(XLON.wrf.T17.raw), as.vector(XLAT.wrf.T17.raw))
colnames(model.T17.df)                        <- c("visib.m", "long", "lat")
#head(model.T17.df)
model.T18.df                                  <- data.frame(as.vector(VISIB.wrf.T18.raw), as.vector(XLON.wrf.T18.raw), as.vector(XLAT.wrf.T18.raw))
colnames(model.T18.df)                        <- c("visib.m", "long", "lat")
#head(model.T18.df)
model.T19.df                                  <- data.frame(as.vector(VISIB.wrf.T19.raw), as.vector(XLON.wrf.T19.raw), as.vector(XLAT.wrf.T19.raw))
colnames(model.T19.df)                        <- c("visib.m", "long", "lat")
#head(model.T19.df)
model.T20.df                                  <- data.frame(as.vector(VISIB.wrf.T20.raw), as.vector(XLON.wrf.T20.raw), as.vector(XLAT.wrf.T20.raw))
colnames(model.T20.df)                        <- c("visib.m", "long", "lat")
#head(model.T20.df)

# build data arrays that make up metar.df
#   NOTE: max allowable visib.m = 9.26 km, to match max allowable model.T#.df
#visib.km           <- c(         9.3,          3.2,          8.1,          4.8,          0.8,          0.8,         1.6,           2.8,          6.5,          6.5,         4.8,           1.6,          1.2,          1.2,          0.8,          1.2,          2.4,          9.3,          9.3,          4.0,          4.8,          8.0,          8.0,          6.5,          6.5,          2.8,          6.5,          3.2,          2.0,          0.4,          0.8,          2.0,          0.4,          2.0,          3.2,          0.8,          1.2,          2.0,          3.2,          2.8,          3.2,          4.8,          6.5,          4.0,          8.1,          9.3,          3.2,          9.3,          2.4,          6.5,          4.0,          4.8,          4.8,          8.0,          9.3,          9.3,          9.3,          3.2,          6.5,          9.3,          9.3,          2.4,          3.2,          9.3,          9.3,          9.3,          9.3,          2.4,          1.6,          4.8,          3.2,          3.2,          4.8,          3.2,          3.2,          3.2,          2.4,          1.6,          1.6,          4.8,          4.0,          3.2,          6.5,          9.3,          9.3,          0.8,          0.4,          0.4,          0.8,          1.2,          2.8,          9.3,          4.0,         9.3,          9.3,           2.8,          3.2,          1.2,          1.6,          1.2,          2.4,          4.8,          1.6,          4.8,          4.8,          1.6,          0.8,          0.4,          0.4,          1.2,          0.8,          1.2,          1.2,          1.2,          8.1,          9.3,          8.1,          9.3,          6.5,          4.0,          4.0,          4.8,          4.8,          4.8,          4.8,          2.0,          1.6,          1.6,          0.8,          3.2,          1.6,          0.4,          0.8,          6.5,          4.8,          8.1,          6.5,          4.0,          6.5,          6.5,          4.8,          3.2,          6.5,          2.8,          6.5,          8.1,          8.1,          3.2,          3.2,          1.6,          4.8,          1.6,          0.8,          1.2,          0.8,          0.8,          1.6,          3.2,          9.3,          3.2,          0.4,          1.6,          9.3,          9.3,          8.1,          2.0,          2.4,          6.5,          9.3,          9.3,          9.3)
#time.HHMM          <- c(     "03:51",      "04:51",      "05:51",      "06:51",      "07:51",      "08:51",     "09:53",       "10:51",      "11:51",      "12:51",     "14:07",       "15:08",      "15:51",      "16:51",      "18:07",      "18:51",      "19:51",      "20:58",      "21:53",      "03:53",      "04:53",      "05:53",      "06:53",      "07:53",      "08:53",      "09:53",      "10:53",      "11:53",      "12:58",      "14:04",      "14:53",      "15:53",      "17:00",      "17:53",      "18:59",      "19:53",      "20:53",      "21:53",      "03:52",      "04:52",      "05:52",      "06:59",      "07:52",      "08:58",      "09:59",      "11:01",      "11:52",      "12:52",      "13:59",      "14:52",      "16:08",      "17:07",      "17:52",      "19:02",      "19:52",      "20:52",      "21:52",      "03:54",      "04:54",      "05:54",      "06:54",      "07:54",      "08:54",      "09:54",      "10:54",      "11:54",      "12:54",      "13:54",      "14:54",      "15:54",      "16:54",      "18:09",      "19:02",      "19:54",      "20:54",      "21:57",      "03:55",      "04:55",      "05:55",      "06:55",      "07:55",      "08:55",      "09:55",      "10:55",      "11:55",      "12:55",      "13:55",      "14:55",      "15:55",      "16:55",      "17:55",      "18:55",      "19:55",     "20:55",      "21:55",       "03:55",      "04:55",      "05:55",      "06:55",      "07:55",      "08:55",      "09:55",      "10:55",      "11:55",      "12:55",      "13:55",      "15:04",      "16:02",      "16:55",      "17:55",      "18:55",      "19:55",      "20:55",      "21:55",      "03:53",      "04:53",      "05:53",      "06:53",      "07:53",      "08:53",      "09:53",      "10:53",      "11:53",      "12:53",      "13:53",      "14:53",      "15:53",      "16:53",      "17:53",      "18:53",      "19:53",      "20:53",      "21:53",      "03:55",      "04:56",      "05:56",      "06:56",      "07:56",      "08:56",      "09:56",      "10:55",      "11:55",      "12:55",      "13:55",      "14:56",      "15:56",      "16:56",      "17:55",      "18:55",      "19:55",      "20:55",      "21:56",      "03:55",      "04:55",      "05:55",      "06:55",      "07:55",      "08:55",      "09:55",      "10:55",      "11:55",      "12:55",      "13:55",      "14:55",      "15:55",      "16:55",      "17:55",      "18:55",      "19:55",      "20:55",      "21:55")
#group              <- c(           4,            5,            6,            7,            8,            9,          10,            11,           12,           13,          14,            15,           16,           17,           18,           19,           20,           21,           22,            4,            5,            6,            7,            8,            9,           10,           11,           12,           13,           14,           15,           16,           17,           18,           19,           20,           21,           22,            4,            5,            6,            7,            8,            9,           10,           11,           12,           13,           14,           15,           16,           17,           18,           19,           20,           21,           22,            4,            5,            6,            7,            8,            9,           10,           11,           12,           13,           14,           15,           16,           17,           18,           19,           20,           21,           22,            4,            5,            6,            7,            8,            9,           10,           11,           12,           13,           14,           15,           16,           17,           18,           19,           20,          21,           22,             4,            5,            6,            7,            8,            9,           10,           11,           12,           13,           14,           15,           16,           17,           18,           19,           20,           21,           22,            4,            5,            6,            7,            8,            9,           10,           11,           12,           13,           14,           15,           16,           17,           18,           19,           20,           21,           22,            4,            5,            6,            7,            8,            9,           10,           11,           12,           13,           14,           15,           16,           17,           18,           19,           20,           21,           22,            4,            5,            6,            7,            8,            9,           10,           11,           12,           13,           14,           15,           16,           17,           18,           19,           20,           21,           22)
#radar              <- c(      "KORD",       "KORD",       "KORD",       "KORD",       "KORD",       "KORD",      "KORD",        "KORD",       "KORD",       "KORD",       "KORD",       "KORD",       "KORD",       "KORD",       "KORD",       "KORD",       "KORD",       "KORD",       "KORD",       "KMWC",       "KMWC",       "KMWC",       "KMWC",       "KMWC",       "KMWC",       "KMWC",       "KMWC",       "KMWC",       "KMWC",       "KMWC",       "KMWC",       "KMWC",       "KMWC",       "KMWC",       "KMWC",       "KMWC",       "KMWC",       "KMWC",       "KDVN",       "KDVN",       "KDVN",       "KDVN",       "KDVN",       "KDVN",       "KDVN",       "KDVN",       "KDVN",       "KDVN",       "KDVN",       "KDVN",       "KDVN",       "KDVN",       "KDVN",       "KDVN",       "KDVN",       "KDVN",       "KDVN",       "KDSM",       "KDSM",       "KDSM",       "KDSM",       "KDSM",       "KDSM",       "KDSM",       "KDSM",       "KDSM",       "KDSM",       "KDSM",       "KDSM",       "KDSM",       "KDSM",       "KDSM",       "KDSM",       "KDSM",       "KDSM",       "KDSM",       "KAAA",       "KAAA",       "KAAA",       "KAAA",       "KAAA",       "KAAA",       "KAAA",       "KAAA",       "KAAA",       "KAAA",       "KAAA",       "KAAA",       "KAAA",       "KAAA",       "KAAA",       "KAAA",       "KAAA",      "KAAA",       "KAAA",        "KDLH",       "KDLH",       "KDLH",       "KDLH",       "KDLH",       "KDLH",       "KDLH",       "KDLH",       "KDLH",       "KDLH",       "KDLH",       "KDLH",       "KDLH",       "KDLH",       "KDLH",       "KDLH",       "KDLH",       "KDLH",       "KDLH",       "KAZO",       "KAZO",       "KAZO",       "KAZO",       "KAZO",       "KAZO",       "KAZO",       "KAZO",       "KAZO",       "KAZO",       "KAZO",       "KAZO",       "KAZO",       "KAZO",       "KAZO",       "KAZO",       "KAZO",       "KAZO",       "KAZO",       "KLWA",       "KLWA",       "KLWA",       "KLWA",       "KLWA",       "KLWA",       "KLWA",       "KLWA",       "KLWA",       "KLWA",       "KLWA",       "KLWA",       "KLWA",       "KLWA",       "KLWA",       "KLWA",       "KLWA",       "KLWA",       "KLWA",       "K3LF",       "K3LF",       "K3LF",       "K3LF",       "K3LF",       "K3LF",       "K3LF",       "K3LF",       "K3LF",       "K3LF",       "K3LF",       "K3LF",       "K3LF",       "K3LF",       "K3LF",       "K3LF",       "K3LF",       "K3LF",       "K3LF")
#lon.deg            <- c(lon.KORD.deg, lon.KORD.deg, lon.KORD.deg, lon.KORD.deg, lon.KORD.deg, lon.KORD.deg, lon.KORD.deg, lon.KORD.deg, lon.KORD.deg, lon.KORD.deg, lon.KORD.deg, lon.KORD.deg, lon.KORD.deg, lon.KORD.deg, lon.KORD.deg, lon.KORD.deg, lon.KORD.deg, lon.KORD.deg, lon.KORD.deg, lon.KMWC.deg, lon.KMWC.deg, lon.KMWC.deg, lon.KMWC.deg, lon.KMWC.deg, lon.KMWC.deg, lon.KMWC.deg, lon.KMWC.deg, lon.KMWC.deg, lon.KMWC.deg, lon.KMWC.deg, lon.KMWC.deg, lon.KMWC.deg, lon.KMWC.deg, lon.KMWC.deg, lon.KMWC.deg, lon.KMWC.deg, lon.KMWC.deg, lon.KMWC.deg, lon.KDVN.deg, lon.KDVN.deg, lon.KDVN.deg, lon.KDVN.deg, lon.KDVN.deg, lon.KDVN.deg, lon.KDVN.deg, lon.KDVN.deg, lon.KDVN.deg, lon.KDVN.deg, lon.KDVN.deg, lon.KDVN.deg, lon.KDVN.deg, lon.KDVN.deg, lon.KDVN.deg, lon.KDVN.deg, lon.KDVN.deg, lon.KDVN.deg, lon.KDVN.deg, lon.KDSM.deg, lon.KDSM.deg, lon.KDSM.deg, lon.KDSM.deg, lon.KDSM.deg, lon.KDSM.deg, lon.KDSM.deg, lon.KDSM.deg, lon.KDSM.deg, lon.KDSM.deg, lon.KDSM.deg, lon.KDSM.deg, lon.KDSM.deg, lon.KDSM.deg, lon.KDSM.deg, lon.KDSM.deg, lon.KDSM.deg, lon.KDSM.deg, lon.KDSM.deg, lon.KAAA.deg, lon.KAAA.deg, lon.KAAA.deg, lon.KAAA.deg, lon.KAAA.deg, lon.KAAA.deg, lon.KAAA.deg, lon.KAAA.deg, lon.KAAA.deg, lon.KAAA.deg, lon.KAAA.deg, lon.KAAA.deg, lon.KAAA.deg, lon.KAAA.deg, lon.KAAA.deg, lon.KAAA.deg, lon.KAAA.deg, lon.KAAA.deg, lon.KAAA.deg, lon.KDLH.deg, lon.KDLH.deg, lon.KDLH.deg, lon.KDLH.deg, lon.KDLH.deg, lon.KDLH.deg, lon.KDLH.deg, lon.KDLH.deg, lon.KDLH.deg, lon.KDLH.deg, lon.KDLH.deg, lon.KDLH.deg, lon.KDLH.deg, lon.KDLH.deg, lon.KDLH.deg, lon.KDLH.deg, lon.KDLH.deg, lon.KDLH.deg, lon.KDLH.deg, lon.KAZO.deg, lon.KAZO.deg, lon.KAZO.deg, lon.KAZO.deg, lon.KAZO.deg, lon.KAZO.deg, lon.KAZO.deg, lon.KAZO.deg, lon.KAZO.deg, lon.KAZO.deg, lon.KAZO.deg, lon.KAZO.deg, lon.KAZO.deg, lon.KAZO.deg, lon.KAZO.deg, lon.KAZO.deg, lon.KAZO.deg, lon.KAZO.deg, lon.KAZO.deg, lon.KLWA.deg, lon.KLWA.deg, lon.KLWA.deg, lon.KLWA.deg, lon.KLWA.deg, lon.KLWA.deg, lon.KLWA.deg, lon.KLWA.deg, lon.KLWA.deg, lon.KLWA.deg, lon.KLWA.deg, lon.KLWA.deg, lon.KLWA.deg, lon.KLWA.deg, lon.KLWA.deg, lon.KLWA.deg, lon.KLWA.deg, lon.KLWA.deg, lon.KLWA.deg, lon.K3LF.deg, lon.K3LF.deg, lon.K3LF.deg, lon.K3LF.deg, lon.K3LF.deg, lon.K3LF.deg, lon.K3LF.deg, lon.K3LF.deg, lon.K3LF.deg, lon.K3LF.deg, lon.K3LF.deg, lon.K3LF.deg, lon.K3LF.deg, lon.K3LF.deg, lon.K3LF.deg, lon.K3LF.deg, lon.K3LF.deg, lon.K3LF.deg, lon.K3LF.deg)
#lat.deg            <- c(lat.KORD.deg, lat.KORD.deg, lat.KORD.deg, lat.KORD.deg, lat.KORD.deg, lat.KORD.deg, lat.KORD.deg, lat.KORD.deg, lat.KORD.deg, lat.KORD.deg, lat.KORD.deg, lat.KORD.deg, lat.KORD.deg, lat.KORD.deg, lat.KORD.deg, lat.KORD.deg, lat.KORD.deg, lat.KORD.deg, lat.KORD.deg, lat.KMWC.deg, lat.KMWC.deg, lat.KMWC.deg, lat.KMWC.deg, lat.KMWC.deg, lat.KMWC.deg, lat.KMWC.deg, lat.KMWC.deg, lat.KMWC.deg, lat.KMWC.deg, lat.KMWC.deg, lat.KMWC.deg, lat.KMWC.deg, lat.KMWC.deg, lat.KMWC.deg, lat.KMWC.deg, lat.KMWC.deg, lat.KMWC.deg, lat.KMWC.deg, lat.KDVN.deg, lat.KDVN.deg, lat.KDVN.deg, lat.KDVN.deg, lat.KDVN.deg, lat.KDVN.deg, lat.KDVN.deg, lat.KDVN.deg, lat.KDVN.deg, lat.KDVN.deg, lat.KDVN.deg, lat.KDVN.deg, lat.KDVN.deg, lat.KDVN.deg, lat.KDVN.deg, lat.KDVN.deg, lat.KDVN.deg, lat.KDVN.deg, lat.KDVN.deg, lat.KDSM.deg, lat.KDSM.deg, lat.KDSM.deg, lat.KDSM.deg, lat.KDSM.deg, lat.KDSM.deg, lat.KDSM.deg, lat.KDSM.deg, lat.KDSM.deg, lat.KDSM.deg, lat.KDSM.deg, lat.KDSM.deg, lat.KDSM.deg, lat.KDSM.deg, lat.KDSM.deg, lat.KDSM.deg, lat.KDSM.deg, lat.KDSM.deg, lat.KDSM.deg, lat.KAAA.deg, lat.KAAA.deg, lat.KAAA.deg, lat.KAAA.deg, lat.KAAA.deg, lat.KAAA.deg, lat.KAAA.deg, lat.KAAA.deg, lat.KAAA.deg, lat.KAAA.deg, lat.KAAA.deg, lat.KAAA.deg, lat.KAAA.deg, lat.KAAA.deg, lat.KAAA.deg, lat.KAAA.deg, lat.KAAA.deg, lat.KAAA.deg, lat.KAAA.deg, lat.KDLH.deg, lat.KDLH.deg, lat.KDLH.deg, lat.KDLH.deg, lat.KDLH.deg, lat.KDLH.deg, lat.KDLH.deg, lat.KDLH.deg, lat.KDLH.deg, lat.KDLH.deg, lat.KDLH.deg, lat.KDLH.deg, lat.KDLH.deg, lat.KDLH.deg, lat.KDLH.deg, lat.KDLH.deg, lat.KDLH.deg, lat.KDLH.deg, lat.KDLH.deg, lat.KAZO.deg, lat.KAZO.deg, lat.KAZO.deg, lat.KAZO.deg, lat.KAZO.deg, lat.KAZO.deg, lat.KAZO.deg, lat.KAZO.deg, lat.KAZO.deg, lat.KAZO.deg, lat.KAZO.deg, lat.KAZO.deg, lat.KAZO.deg, lat.KAZO.deg, lat.KAZO.deg, lat.KAZO.deg, lat.KAZO.deg, lat.KAZO.deg, lat.KAZO.deg, lat.KLWA.deg, lat.KLWA.deg, lat.KLWA.deg, lat.KLWA.deg, lat.KLWA.deg, lat.KLWA.deg, lat.KLWA.deg, lat.KLWA.deg, lat.KLWA.deg, lat.KLWA.deg, lat.KLWA.deg, lat.KLWA.deg, lat.KLWA.deg, lat.KLWA.deg, lat.KLWA.deg, lat.KLWA.deg, lat.KLWA.deg, lat.KLWA.deg, lat.KLWA.deg, lat.K3LF.deg, lat.K3LF.deg, lat.K3LF.deg, lat.K3LF.deg, lat.K3LF.deg, lat.K3LF.deg, lat.K3LF.deg, lat.K3LF.deg, lat.K3LF.deg, lat.K3LF.deg, lat.K3LF.deg, lat.K3LF.deg, lat.K3LF.deg, lat.K3LF.deg, lat.K3LF.deg, lat.K3LF.deg, lat.K3LF.deg, lat.K3LF.deg, lat.K3LF.deg)
#visib.m            <- visib.km * 1000
visib.T2.km           <- c(         9.3,           4.0,          3.2,          3.2,          2.4,          2.8,          8.1,          6.5,          0.8)
radar.T2              <- c(      "KORD",        "KMWC",       "KDVN",       "KDSM",       "KAAA",       "KDLH",       "KAZO",       "KLWA",       "K3LF")
lon.T2.deg            <- c(lon.KORD.deg,  lon.KMWC.deg, lon.KDVN.deg, lon.KDSM.deg, lon.KAAA.deg, lon.KDLH.deg, lon.KAZO.deg, lon.KLWA.deg, lon.K3LF.deg)
lat.T2.deg            <- c(lat.KORD.deg,  lat.KMWC.deg, lat.KDVN.deg, lat.KDSM.deg, lat.KAAA.deg, lat.KDLH.deg, lat.KAZO.deg, lat.KLWA.deg, lat.K3LF.deg)
visib.T2.m            <- visib.T2.km * 1000
visib.T3.km           <- c(         3.2,           4.8,          2.8,          6.5,          1.6,          3.2,          9.3,          4.8,          1.2)
radar.T3              <- c(      "KORD",        "KMWC",       "KDVN",       "KDSM",       "KAAA",       "KDLH",       "KAZO",       "KLWA",       "K3LF")
lon.T3.deg            <- c(lon.KORD.deg,  lon.KMWC.deg, lon.KDVN.deg, lon.KDSM.deg, lon.KAAA.deg, lon.KDLH.deg, lon.KAZO.deg, lon.KLWA.deg, lon.K3LF.deg)
lat.T3.deg            <- c(lat.KORD.deg,  lat.KMWC.deg, lat.KDVN.deg, lat.KDSM.deg, lat.KAAA.deg, lat.KDLH.deg, lat.KAZO.deg, lat.KLWA.deg, lat.K3LF.deg)
visib.T3.m            <- visib.T3.km * 1000
visib.T4.km           <- c(         8.1,           8.1,          3.2,          9.3,          1.6,          1.2,          8.1,          8.1,          0.8)
radar.T4              <- c(      "KORD",        "KMWC",       "KDVN",       "KDSM",       "KAAA",       "KDLH",       "KAZO",       "KLWA",       "K3LF")
lon.T4.deg            <- c(lon.KORD.deg,  lon.KMWC.deg, lon.KDVN.deg, lon.KDSM.deg, lon.KAAA.deg, lon.KDLH.deg, lon.KAZO.deg, lon.KLWA.deg, lon.K3LF.deg)
lat.T4.deg            <- c(lat.KORD.deg,  lat.KMWC.deg, lat.KDVN.deg, lat.KDSM.deg, lat.KAAA.deg, lat.KDLH.deg, lat.KAZO.deg, lat.KLWA.deg, lat.K3LF.deg)
visib.T4.m            <- visib.T4.km * 1000
visib.T5.km           <- c(         4.8,           8.1,          4.8,          9.3,          4.8,          1.6,          9.3,          6.5,          0.8)
radar.T5              <- c(      "KORD",        "KMWC",       "KDVN",       "KDSM",       "KAAA",       "KDLH",       "KAZO",       "KLWA",       "K3LF")
lon.T5.deg            <- c(lon.KORD.deg,  lon.KMWC.deg, lon.KDVN.deg, lon.KDSM.deg, lon.KAAA.deg, lon.KDLH.deg, lon.KAZO.deg, lon.KLWA.deg, lon.K3LF.deg)
lat.T5.deg            <- c(lat.KORD.deg,  lat.KMWC.deg, lat.KDVN.deg, lat.KDSM.deg, lat.KAAA.deg, lat.KDLH.deg, lat.KAZO.deg, lat.KLWA.deg, lat.K3LF.deg)
visib.T5.m            <- visib.T5.km * 1000
visib.T6.km           <- c(         0.8,           6.5,          6.5,          2.4,          4.0,          1.2,          6.5,          4.0,          1.6)
radar.T6              <- c(      "KORD",        "KMWC",       "KDVN",       "KDSM",       "KAAA",       "KDLH",       "KAZO",       "KLWA",       "K3LF")
lon.T6.deg            <- c(lon.KORD.deg,  lon.KMWC.deg, lon.KDVN.deg, lon.KDSM.deg, lon.KAAA.deg, lon.KDLH.deg, lon.KAZO.deg, lon.KLWA.deg, lon.K3LF.deg)
lat.T6.deg            <- c(lat.KORD.deg,  lat.KMWC.deg, lat.KDVN.deg, lat.KDSM.deg, lat.KAAA.deg, lat.KDLH.deg, lat.KAZO.deg, lat.KLWA.deg, lat.K3LF.deg)
visib.T6.m            <- visib.T6.km * 1000
visib.T7.km           <- c(         0.8,           6.5,          4.0,          3.2,          3.2,          2.4,          4.0,          6.5,          3.2)
radar.T7              <- c(      "KORD",        "KMWC",       "KDVN",       "KDSM",       "KAAA",       "KDLH",       "KAZO",       "KLWA",       "K3LF")
lon.T7.deg            <- c(lon.KORD.deg,  lon.KMWC.deg, lon.KDVN.deg, lon.KDSM.deg, lon.KAAA.deg, lon.KDLH.deg, lon.KAZO.deg, lon.KLWA.deg, lon.K3LF.deg)
lat.T7.deg            <- c(lat.KORD.deg,  lat.KMWC.deg, lat.KDVN.deg, lat.KDSM.deg, lat.KAAA.deg, lat.KDLH.deg, lat.KAZO.deg, lat.KLWA.deg, lat.K3LF.deg)
visib.T7.m            <- visib.T7.km * 1000
visib.T8.km           <- c(         1.6,           2.8,          8.1,          9.3,          6.5,          4.8,          4.0,          6.5,          9.3)
radar.T8              <- c(      "KORD",        "KMWC",       "KDVN",       "KDSM",       "KAAA",       "KDLH",       "KAZO",       "KLWA",       "K3LF")
lon.T8.deg            <- c(lon.KORD.deg,  lon.KMWC.deg, lon.KDVN.deg, lon.KDSM.deg, lon.KAAA.deg, lon.KDLH.deg, lon.KAZO.deg, lon.KLWA.deg, lon.K3LF.deg)
lat.T8.deg            <- c(lat.KORD.deg,  lat.KMWC.deg, lat.KDVN.deg, lat.KDSM.deg, lat.KAAA.deg, lat.KDLH.deg, lat.KAZO.deg, lat.KLWA.deg, lat.K3LF.deg)
visib.T8.m            <- visib.T8.km * 1000
visib.T9.km           <- c(         2.8,           6.5,          9.3,          9.3,          9.3,          1.6,          4.8,          4.8,          3.2)
radar.T9              <- c(      "KORD",        "KMWC",       "KDVN",       "KDSM",       "KAAA",       "KDLH",       "KAZO",       "KLWA",       "K3LF")
lon.T9.deg            <- c(lon.KORD.deg,  lon.KMWC.deg, lon.KDVN.deg, lon.KDSM.deg, lon.KAAA.deg, lon.KDLH.deg, lon.KAZO.deg, lon.KLWA.deg, lon.K3LF.deg)
lat.T9.deg            <- c(lat.KORD.deg,  lat.KMWC.deg, lat.KDVN.deg, lat.KDSM.deg, lat.KAAA.deg, lat.KDLH.deg, lat.KAZO.deg, lat.KLWA.deg, lat.K3LF.deg)
visib.T9.m            <- visib.T9.km * 1000
visib.T10.km           <- c(         6.5,           3.2,          9.3,          9.3,          9.3,          4.8,          2.4,          3.2,          0.4)
radar.T10              <- c(      "KORD",        "KMWC",       "KDVN",       "KDSM",       "KAAA",       "KDLH",       "KAZO",       "KLWA",       "K3LF")
lon.T10.deg            <- c(lon.KORD.deg,  lon.KMWC.deg, lon.KDVN.deg, lon.KDSM.deg, lon.KAAA.deg, lon.KDLH.deg, lon.KAZO.deg, lon.KLWA.deg, lon.K3LF.deg)
lat.T10.deg            <- c(lat.KORD.deg,  lat.KMWC.deg, lat.KDVN.deg, lat.KDSM.deg, lat.KAAA.deg, lat.KDLH.deg, lat.KAZO.deg, lat.KLWA.deg, lat.K3LF.deg)
visib.T10.m            <- visib.T10.km * 1000
visib.T11.km           <- c(         6.5,           2.0,          9.3,          9.3,          0.8,          4.0,          4.8,          6.5,          1.6)
radar.T11              <- c(      "KORD",        "KMWC",       "KDVN",       "KDSM",       "KAAA",       "KDLH",       "KAZO",       "KLWA",       "K3LF")
lon.T11.deg            <- c(lon.KORD.deg,  lon.KMWC.deg, lon.KDVN.deg, lon.KDSM.deg, lon.KAAA.deg, lon.KDLH.deg, lon.KAZO.deg, lon.KLWA.deg, lon.K3LF.deg)
lat.T11.deg            <- c(lat.KORD.deg,  lat.KMWC.deg, lat.KDVN.deg, lat.KDSM.deg, lat.KAAA.deg, lat.KDLH.deg, lat.KAZO.deg, lat.KLWA.deg, lat.K3LF.deg)
visib.T11.m            <- visib.T11.km * 1000
visib.T12.km           <- c(         4.8,           0.8,          3.2,          2.4,          0.4,          1.6,          4.8,          2.8,          9.3)
radar.T12              <- c(      "KORD",        "KMWC",       "KDVN",       "KDSM",       "KAAA",       "KDLH",       "KAZO",       "KLWA",       "K3LF")
lon.T12.deg            <- c(lon.KORD.deg,  lon.KMWC.deg, lon.KDVN.deg, lon.KDSM.deg, lon.KAAA.deg, lon.KDLH.deg, lon.KAZO.deg, lon.KLWA.deg, lon.K3LF.deg)
lat.T12.deg            <- c(lat.KORD.deg,  lat.KMWC.deg, lat.KDVN.deg, lat.KDSM.deg, lat.KAAA.deg, lat.KDLH.deg, lat.KAZO.deg, lat.KLWA.deg, lat.K3LF.deg)
visib.T12.m            <- visib.T12.km * 1000
visib.T13.km           <- c(         1.6,           0.8,          6.5,          1.6,          0.4,          2.4,          4.8,          6.5,          9.3)
radar.T13              <- c(      "KORD",        "KMWC",       "KDVN",       "KDSM",       "KAAA",       "KDLH",       "KAZO",       "KLWA",       "K3LF")
lon.T13.deg            <- c(lon.KORD.deg,  lon.KMWC.deg, lon.KDVN.deg, lon.KDSM.deg, lon.KAAA.deg, lon.KDLH.deg, lon.KAZO.deg, lon.KLWA.deg, lon.K3LF.deg)
lat.T13.deg            <- c(lat.KORD.deg,  lat.KMWC.deg, lat.KDVN.deg, lat.KDSM.deg, lat.KAAA.deg, lat.KDLH.deg, lat.KAZO.deg, lat.KLWA.deg, lat.K3LF.deg)
visib.T13.m            <- visib.T13.km * 1000
visib.T14.km           <- c(         1.2,           2.0,          4.0,          4.8,          0.8,          0.4,          1.6,          8.1,          8.1)
radar.T14              <- c(      "KORD",        "KMWC",       "KDVN",       "KDSM",       "KAAA",       "KDLH",       "KAZO",       "KLWA",       "K3LF")
lon.T14.deg            <- c(lon.KORD.deg,  lon.KMWC.deg, lon.KDVN.deg, lon.KDSM.deg, lon.KAAA.deg, lon.KDLH.deg, lon.KAZO.deg, lon.KLWA.deg, lon.K3LF.deg)
lat.T14.deg            <- c(lat.KORD.deg,  lat.KMWC.deg, lat.KDVN.deg, lat.KDSM.deg, lat.KAAA.deg, lat.KDLH.deg, lat.KAZO.deg, lat.KLWA.deg, lat.K3LF.deg)
visib.T14.m            <- visib.T14.km * 1000
visib.T15.km           <- c(         1.2,           0.4,          4.0,          3.2,          1.2,          0.4,          1.6,          8.1,          2.0)
radar.T15              <- c(      "KORD",        "KMWC",       "KDVN",       "KDSM",       "KAAA",       "KDLH",       "KAZO",       "KLWA",       "K3LF")
lon.T15.deg            <- c(lon.KORD.deg,  lon.KMWC.deg, lon.KDVN.deg, lon.KDSM.deg, lon.KAAA.deg, lon.KDLH.deg, lon.KAZO.deg, lon.KLWA.deg, lon.K3LF.deg)
lat.T15.deg            <- c(lat.KORD.deg,  lat.KMWC.deg, lat.KDVN.deg, lat.KDSM.deg, lat.KAAA.deg, lat.KDLH.deg, lat.KAZO.deg, lat.KLWA.deg, lat.K3LF.deg)
visib.T15.m            <- visib.T15.km * 1000
visib.T16.km           <- c(         0.8,           2.0,          4.8,          6.5,          2.8,          1.2,          0.8,          3.2,          2.4)
radar.T16              <- c(      "KORD",        "KMWC",       "KDVN",       "KDSM",       "KAAA",       "KDLH",       "KAZO",       "KLWA",       "K3LF")
lon.T16.deg            <- c(lon.KORD.deg,  lon.KMWC.deg, lon.KDVN.deg, lon.KDSM.deg, lon.KAAA.deg, lon.KDLH.deg, lon.KAZO.deg, lon.KLWA.deg, lon.K3LF.deg)
lat.T16.deg            <- c(lat.KORD.deg,  lat.KMWC.deg, lat.KDVN.deg, lat.KDSM.deg, lat.KAAA.deg, lat.KDLH.deg, lat.KAZO.deg, lat.KLWA.deg, lat.K3LF.deg)
visib.T16.m            <- visib.T16.km * 1000
visib.T17.km           <- c(         1.2,           3.2,          8.1,          3.2,          9.3,          0.8,          3.2,          3.2,          6.5)
radar.T17              <- c(      "KORD",        "KMWC",       "KDVN",       "KDSM",       "KAAA",       "KDLH",       "KAZO",       "KLWA",       "K3LF")
lon.T17.deg            <- c(lon.KORD.deg,  lon.KMWC.deg, lon.KDVN.deg, lon.KDSM.deg, lon.KAAA.deg, lon.KDLH.deg, lon.KAZO.deg, lon.KLWA.deg, lon.K3LF.deg)
lat.T17.deg            <- c(lat.KORD.deg,  lat.KMWC.deg, lat.KDVN.deg, lat.KDSM.deg, lat.KAAA.deg, lat.KDLH.deg, lat.KAZO.deg, lat.KLWA.deg, lat.K3LF.deg)
visib.T17.m            <- visib.T17.km * 1000
visib.T18.km           <- c(         2.4,           0.8,          9.3,          3.2,          4.0,          1.2,          0.4,          1.6,          9.3)
radar.T18              <- c(      "KORD",        "KMWC",       "KDVN",       "KDSM",       "KAAA",       "KDLH",       "KAZO",       "KLWA",       "K3LF")
lon.T18.deg            <- c(lon.KORD.deg,  lon.KMWC.deg, lon.KDVN.deg, lon.KDSM.deg, lon.KAAA.deg, lon.KDLH.deg, lon.KAZO.deg, lon.KLWA.deg, lon.K3LF.deg)
lat.T18.deg            <- c(lat.KORD.deg,  lat.KMWC.deg, lat.KDVN.deg, lat.KDSM.deg, lat.KAAA.deg, lat.KDLH.deg, lat.KAZO.deg, lat.KLWA.deg, lat.K3LF.deg)
visib.T18.m            <- visib.T18.km * 1000
visib.T19.km           <- c(         9.3,           1.2,          9.3,          3.2,          9.3,          1.2,          0.4,          4.8,          9.3)
radar.T19              <- c(      "KORD",        "KMWC",       "KDVN",       "KDSM",       "KAAA",       "KDLH",       "KAZO",       "KLWA",       "K3LF")
lon.T19.deg            <- c(lon.KORD.deg,  lon.KMWC.deg, lon.KDVN.deg, lon.KDSM.deg, lon.KAAA.deg, lon.KDLH.deg, lon.KAZO.deg, lon.KLWA.deg, lon.K3LF.deg)
lat.T19.deg            <- c(lat.KORD.deg,  lat.KMWC.deg, lat.KDVN.deg, lat.KDSM.deg, lat.KAAA.deg, lat.KDLH.deg, lat.KAZO.deg, lat.KLWA.deg, lat.K3LF.deg)
visib.T19.m            <- visib.T19.km * 1000
visib.T20.km           <- c(         9.3,           2.0,          9.3,          3.2,          9.3,          1.2,          0.8,          1.6,          9.3)
radar.T20              <- c(      "KORD",        "KMWC",       "KDVN",       "KDSM",       "KAAA",       "KDLH",       "KAZO",       "KLWA",       "K3LF")
lon.T20.deg            <- c(lon.KORD.deg,  lon.KMWC.deg, lon.KDVN.deg, lon.KDSM.deg, lon.KAAA.deg, lon.KDLH.deg, lon.KAZO.deg, lon.KLWA.deg, lon.K3LF.deg)
lat.T20.deg            <- c(lat.KORD.deg,  lat.KMWC.deg, lat.KDVN.deg, lat.KDSM.deg, lat.KAAA.deg, lat.KDLH.deg, lat.KAZO.deg, lat.KLWA.deg, lat.K3LF.deg)
visib.T20.m            <- visib.T20.km * 1000

# build metar df, over 9 locations x 20 hours
metar.T2.df           <- data.frame(radar.T2, lon.T2.deg, lat.T2.deg, visib.T2.m)
colnames(metar.T2.df) <- c("radar", "lon.deg", "lat.deg", "visib.m")
head(metar.T2.df)
metar.T3.df           <- data.frame(radar.T3, lon.T3.deg, lat.T3.deg, visib.T3.m)
colnames(metar.T3.df) <- c("radar", "lon.deg", "lat.deg", "visib.m")
head(metar.T3.df)
metar.T4.df           <- data.frame(radar.T4, lon.T4.deg, lat.T4.deg, visib.T4.m)
colnames(metar.T4.df) <- c("radar", "lon.deg", "lat.deg", "visib.m")
head(metar.T4.df)
metar.T5.df           <- data.frame(radar.T5, lon.T5.deg, lat.T5.deg, visib.T5.m)
colnames(metar.T5.df) <- c("radar", "lon.deg", "lat.deg", "visib.m")
head(metar.T5.df)
metar.T6.df           <- data.frame(radar.T6, lon.T6.deg, lat.T6.deg, visib.T6.m)
colnames(metar.T6.df) <- c("radar", "lon.deg", "lat.deg", "visib.m")
head(metar.T6.df)
metar.T7.df           <- data.frame(radar.T7, lon.T7.deg, lat.T7.deg, visib.T7.m)
colnames(metar.T7.df) <- c("radar", "lon.deg", "lat.deg", "visib.m")
head(metar.T7.df)
metar.T8.df           <- data.frame(radar.T8, lon.T8.deg, lat.T8.deg, visib.T8.m)
colnames(metar.T8.df) <- c("radar", "lon.deg", "lat.deg", "visib.m")
head(metar.T8.df)
metar.T9.df           <- data.frame(radar.T9, lon.T9.deg, lat.T9.deg, visib.T9.m)
colnames(metar.T9.df) <- c("radar", "lon.deg", "lat.deg", "visib.m")
head(metar.T9.df)
metar.T10.df           <- data.frame(radar.T10, lon.T10.deg, lat.T10.deg, visib.T10.m)
colnames(metar.T10.df) <- c("radar", "lon.deg", "lat.deg", "visib.m")
head(metar.T10.df)
metar.T11.df           <- data.frame(radar.T11, lon.T11.deg, lat.T11.deg, visib.T11.m)
colnames(metar.T11.df) <- c("radar", "lon.deg", "lat.deg", "visib.m")
head(metar.T11.df)
metar.T12.df           <- data.frame(radar.T12, lon.T12.deg, lat.T12.deg, visib.T12.m)
colnames(metar.T12.df) <- c("radar", "lon.deg", "lat.deg", "visib.m")
head(metar.T12.df)
metar.T13.df           <- data.frame(radar.T13, lon.T13.deg, lat.T13.deg, visib.T13.m)
colnames(metar.T13.df) <- c("radar", "lon.deg", "lat.deg", "visib.m")
head(metar.T13.df)
metar.T14.df           <- data.frame(radar.T14, lon.T14.deg, lat.T14.deg, visib.T14.m)
colnames(metar.T14.df) <- c("radar", "lon.deg", "lat.deg", "visib.m")
head(metar.T14.df)
metar.T15.df           <- data.frame(radar.T15, lon.T15.deg, lat.T15.deg, visib.T15.m)
colnames(metar.T15.df) <- c("radar", "lon.deg", "lat.deg", "visib.m")
head(metar.T15.df)
metar.T16.df           <- data.frame(radar.T16, lon.T16.deg, lat.T16.deg, visib.T16.m)
colnames(metar.T16.df) <- c("radar", "lon.deg", "lat.deg", "visib.m")
head(metar.T16.df)
metar.T17.df           <- data.frame(radar.T17, lon.T17.deg, lat.T17.deg, visib.T17.m)
colnames(metar.T17.df) <- c("radar", "lon.deg", "lat.deg", "visib.m")
head(metar.T17.df)
metar.T18.df           <- data.frame(radar.T18, lon.T18.deg, lat.T18.deg, visib.T18.m)
colnames(metar.T18.df) <- c("radar", "lon.deg", "lat.deg", "visib.m")
head(metar.T18.df)
metar.T19.df           <- data.frame(radar.T19, lon.T19.deg, lat.T19.deg, visib.T19.m)
colnames(metar.T19.df) <- c("radar", "lon.deg", "lat.deg", "visib.m")
head(metar.T19.df)
metar.T20.df           <- data.frame(radar.T20, lon.T20.deg, lat.T20.deg, visib.T20.m)
colnames(metar.T20.df) <- c("radar", "lon.deg", "lat.deg", "visib.m")
head(metar.T20.df)

#-----------------------------
# compute some statistics of model vis relative to metar vis
#-----------------------------
diff.metar.v.model.vis.m        <- zeros(dim(data.metar.visib)[1], 1)
ind.closest.time                <- zeros(dim(data.metar.visib)[1], 1)
for (i in 1: length(data.metar.visib$HH.MM)) {
  ind.closest.time[i]         <- which.min(abs(data.metar.visib$HH.MM[i] - data.model.visib$HH.MM))
  diff.metar.v.model.vis.m[i] <- data.metar.visib$visib.m[i] - data.model.visib$visib.m[ind.closest.time]
}
#mean.diff.metar.v.model.vis.m   <- mean(diff.metar.v.model.vis.m)
#median.diff.metar.v.model.vis.m <- median(diff.metar.v.model.vis.m)
#std.diff.metar.v.model.vis.m    <- std(diff.metar.v.model.vis.m)

#print(paste("The mean diff in model v metar vis was ", round(mean.diff.metar.v.model.vis.m, digits=0), " m", sep=""))
#print(paste("The median diff in model v metar vis was ", round(median.diff.metar.v.model.vis.m, digits=0), " m", sep=""))
#print(paste("The std diff in model v metar vis was ", round(std.diff.metar.v.model.vis.m, digits=0), " m", sep=""))

#-----------------------------
# plotting
#-----------------------------
# plot vis-H derived from model versus vis-H and RVR from METARs
plot <- ggplot(data=data.combine.visib, aes(x=combine.HH.MM, y=combine.visib, group=combine.obs)) + 
#        geom_line(aes(color=combine.obs)) +
        geom_point(aes(color=combine.obs)) + 
        theme(legend.position="right") +
        labs(title=paste("2019-02-07 VIS @ ", metar.station, sep=""), x="Time [UTC]", y="VIS [m]") +
        ylim(min(data.combine.visib$combine.visib)-50, 50000) 
plot
#plot + scale_y_continuous(trans='log10')

# scatterplot
plot(data.model.visib$visib.m[ind.closest.time], data.metar.visib$visib.m, xlim=c(0,20000), ylim=c(0,20000), xlab="VIS-model [m]", ylab="VIS-METAR [m]", main=paste("2019-02-07, 3:00 to 22:00 UTC @ ", metar.station, " Model-to-METAR VIS Correlation", sep="") )
abline(0,1, lwd=2, col="red")
grid()
#lines(lowess(data.model.visib$visib.m[ind.closest.time], data.metar.visib$visib.m), col="blue")
pearson.coeff   <- round(cor(data.model.visib$visib.m[ind.closest.time], data.metar.visib$visib.m, method="pearson", use="complete.obs"), digits=2)
num.pts.pearson <- length(ind.closest.time)
text(8000, 16000, paste("Pearson correlation = R = ", pearson.coeff, sep=""))
#text(8000, 16000, paste("Pearson correlation = R = ", pearson.coeff, " (mod + correlation)", sep=""))
text(8000, 15000, paste("R^2 = ", round(pearson.coeff^2, digits=2), sep=""))
text(8000, 14000, paste("Num pts = ", num.pts.pearson, sep=""))
legend(2000, 25000, legend=c("1:1 line", "Lowess best-fit"), col=c("red", "blue"), lty=c(1,1))

states                             <- map_data("state")
icicle.domain                      <- subset(states, region %in% c("minnesota", "iowa", "missouri", "illinois", "wisconsin",  "michigan", "indiana", "ohio", "kentucky"))

# ggplot of model Time #2 visib, with midwest map and color filled circles from METAR
ggplot(data = model.T2.df, aes(x=long, y=lat)) +  
  geom_point(aes(color = visib.m)) +
  geom_point(data=icicle.domain, aes(x=long, y=lat), fill=NA, color="grey50", size=0.25) +
  xlim(-98, -85) +
  ylim( 36,  50) +
  geom_point(data=metar.T2.df, aes(x=lon.deg, y=lat.deg, fill=visib.m), shape=21, color="white", size=3, stroke=1) +
  geom_text( data=metar.T2.df, aes(x=lon.deg, y=lat.deg, label=radar),            color="white", hjust=0.5, vjust=2.0) +
  coord_fixed() +
  #scale_color_brewer(palette = "PRGn") +
  ggtitle("2019/02/07 @ 04:00 UTC, WRF 3-km VIS (K-G)") +
  theme(legend.position      = "bottom") +
  theme(legend.key.width     = unit(2, "cm"))

# ggplot of model Time #3 1 visib, with midwest map and color filled circles from METAR
ggplot(data = model.T3.df, aes(x=long, y=lat)) +  
  geom_point(aes(color = visib.m)) +
  geom_point(data=icicle.domain, aes(x=long, y=lat), fill=NA, color="grey50", size=0.25) +
  xlim(-98, -85) +
  ylim( 36,  50) +
  geom_point(data=metar.T3.df, aes(x=lon.deg, y=lat.deg, fill=visib.m), shape=21, color="white", size=3, stroke=1) +
  geom_text( data=metar.T3.df, aes(x=lon.deg, y=lat.deg, label=radar),            color="white", hjust=0.5, vjust=2.0) +
  coord_fixed() +
  #scale_color_brewer(palette = "PRGn") +
  ggtitle("2019/02/07 @ 05:00 UTC, WRF 3-km VIS (K-G)") +
  theme(legend.position      = "bottom") +
  theme(legend.key.width     = unit(2, "cm"))

# ggplot of model Time #4 1 visib, with midwest map and color filled circles from METAR
ggplot(data = model.T4.df, aes(x=long, y=lat)) +  
  geom_point(aes(color = visib.m)) +
  geom_point(data=icicle.domain, aes(x=long, y=lat), fill=NA, color="grey50", size=0.25) +
  xlim(-98, -85) +
  ylim( 36,  50) +
  geom_point(data=metar.T4.df, aes(x=lon.deg, y=lat.deg, fill=visib.m), shape=21, color="white", size=3, stroke=1) +
  geom_text( data=metar.T4.df, aes(x=lon.deg, y=lat.deg, label=radar),            color="white", hjust=0.5, vjust=2.0) +
  coord_fixed() +
  #scale_color_brewer(palette = "PRGn") +
  ggtitle("2019/02/07 @ 06:00 UTC, WRF 3-km VIS (K-G)") +
  theme(legend.position      = "bottom") +
  theme(legend.key.width     = unit(2, "cm"))

# ggplot of model Time #5 1 visib, with midwest map and color filled circles from METAR
ggplot(data = model.T5.df, aes(x=long, y=lat)) +  
  geom_point(aes(color = visib.m)) +
  geom_point(data=icicle.domain, aes(x=long,    y=lat), fill=NA, color="grey50", size=0.25) +
  xlim(-98, -85) +
  ylim( 36,  50) +
  geom_point(data=metar.T5.df,      aes(x=lon.deg, y=lat.deg, fill=visib.m), shape=21, color="white", size=3, stroke=1) +
  geom_text( data=metar.T5.df,      aes(x=lon.deg, y=lat.deg, label=radar),            color="white", hjust=0.5, vjust=2.0) +
  coord_fixed() +
  #scale_color_brewer(palette = "PRGn") +
  ggtitle("2019/02/07 @ 07:00 UTC, WRF 3-km VIS (K-G)") +
  theme(legend.position      = "bottom") +
  theme(legend.key.width     = unit(2, "cm"))

# ggplot of model Time #6 1 visib, with midwest map and color filled circles from METAR
ggplot(data = model.T6.df, aes(x=long, y=lat)) +  
  geom_point(aes(color = visib.m)) +
  geom_point(data=icicle.domain, aes(x=long, y=lat), fill=NA, color="grey50", size=0.25) +
  xlim(-98, -85) +
  ylim( 36,  50) +
  geom_point(data=metar.T6.df, aes(x=lon.deg, y=lat.deg, fill=visib.m), shape=21, color="white", size=3, stroke=1) +
  geom_text( data=metar.T6.df, aes(x=lon.deg, y=lat.deg, label=radar),            color="white", hjust=0.5, vjust=2.0) +
  coord_fixed() +
  #scale_color_brewer(palette = "PRGn") +
  ggtitle("2019/02/07 @ 08:00 UTC, WRF 3-km VIS (K-G)") +
  theme(legend.position      = "bottom") +
  theme(legend.key.width     = unit(2, "cm"))

# ggplot of model Time #7 1 visib, with midwest map and color filled circles from METAR
ggplot(data = model.T7.df, aes(x=long, y=lat)) +  
  geom_point(aes(color = visib.m)) +
  geom_point(data=icicle.domain, aes(x=long, y=lat), fill=NA, color="grey50", size=0.25) +
  xlim(-98, -85) +
  ylim( 36,  50) +
  geom_point(data=metar.T7.df, aes(x=lon.deg, y=lat.deg, fill=visib.m), shape=21, color="white", size=3, stroke=1) +
  geom_text( data=metar.T7.df, aes(x=lon.deg, y=lat.deg, label=radar),            color="white", hjust=0.5, vjust=2.0) +
  coord_fixed() +
  #scale_color_brewer(palette = "PRGn") +
  ggtitle("2019/02/07 @ 09:00 UTC, WRF 3-km VIS (K-G)") +
  theme(legend.position      = "bottom") +
  theme(legend.key.width     = unit(2, "cm"))

# ggplot of model Time #8 1 visib, with midwest map and color filled circles from METAR
ggplot(data = model.T8.df, aes(x=long, y=lat)) +  
  geom_point(aes(color = visib.m)) +
  geom_point(data=icicle.domain, aes(x=long, y=lat), fill=NA, color="grey50", size=0.25) +
  xlim(-98, -85) +
  ylim( 36,  50) +
  geom_point(data=metar.T8.df, aes(x=lon.deg, y=lat.deg, fill=visib.m), shape=21, color="white", size=3, stroke=1) +
  geom_text( data=metar.T8.df, aes(x=lon.deg, y=lat.deg, label=radar),            color="white", hjust=0.5, vjust=2.0) +
  coord_fixed() +
  #scale_color_brewer(palette = "PRGn") +
  ggtitle("2019/02/07 @ 10:00 UTC, WRF 3-km VIS (K-G)") +
  theme(legend.position      = "bottom") +
  theme(legend.key.width     = unit(2, "cm"))

# ggplot of model Time #9 1 visib, with midwest map and color filled circles from METAR
ggplot(data = model.T9.df, aes(x=long, y=lat)) +  
  geom_point(aes(color = visib.m)) +
  geom_point(data=icicle.domain, aes(x=long, y=lat), fill=NA, color="grey50", size=0.25) +
  xlim(-98, -85) +
  ylim( 36,  50) +
  geom_point(data=metar.T9.df, aes(x=lon.deg, y=lat.deg, fill=visib.m), shape=21, color="white", size=3, stroke=1) +
  geom_text( data=metar.T9.df, aes(x=lon.deg, y=lat.deg, label=radar),            color="white", hjust=0.5, vjust=2.0) +
  coord_fixed() +
  #scale_color_brewer(palette = "PRGn") +
  ggtitle("2019/02/07 @ 11:00 UTC, WRF 3-km VIS (K-G)") +
  theme(legend.position      = "bottom") +
  theme(legend.key.width     = unit(2, "cm"))

# ggplot of model Time #10 1 visib, with midwest map and color filled circles from METAR
ggplot(data = model.T10.df, aes(x=long, y=lat)) +  
  geom_point(aes(color = visib.m)) +
  geom_point(data=icicle.domain, aes(x=long, y=lat), fill=NA, color="grey50", size=0.25) +
  xlim(-98, -85) +
  ylim( 36,  50) +
  geom_point(data=metar.T10.df, aes(x=lon.deg, y=lat.deg, fill=visib.m), shape=21, color="white", size=3, stroke=1) +
  geom_text( data=metar.T10.df, aes(x=lon.deg, y=lat.deg, label=radar),            color="white", hjust=0.5, vjust=2.0) +
  coord_fixed() +
  #scale_color_brewer(palette = "PRGn") +
  ggtitle("2019/02/07 @ 12:00 UTC, WRF 3-km VIS (K-G)") +
  theme(legend.position      = "bottom") +
  theme(legend.key.width     = unit(2, "cm"))

# ggplot of model Time #11 1 visib, with midwest map and color filled circles from METAR
ggplot(data = model.T11.df, aes(x=long, y=lat)) +  
  geom_point(aes(color = visib.m)) +
  geom_point(data=icicle.domain, aes(x=long, y=lat), fill=NA, color="grey50", size=0.25) +
  xlim(-98, -85) +
  ylim( 36,  50) +
  geom_point(data=metar.T11.df, aes(x=lon.deg, y=lat.deg, fill=visib.m), shape=21, color="white", size=3, stroke=1) +
  geom_text( data=metar.T11.df, aes(x=lon.deg, y=lat.deg, label=radar),            color="white", hjust=0.5, vjust=2.0) +
  coord_fixed() +
  #scale_color_brewer(palette = "PRGn") +
  ggtitle("2019/02/07 @ 13:00 UTC, WRF 3-km VIS (K-G)") +
  theme(legend.position      = "bottom") +
  theme(legend.key.width     = unit(2, "cm"))

# ggplot of model Time #12 1 visib, with midwest map and color filled circles from METAR
ggplot(data = model.T12.df, aes(x=long, y=lat)) +  
  geom_point(aes(color = visib.m)) +
  geom_point(data=icicle.domain, aes(x=long, y=lat), fill=NA, color="grey50", size=0.25) +
  xlim(-98, -85) +
  ylim( 36,  50) +
  geom_point(data=metar.T12.df, aes(x=lon.deg, y=lat.deg, fill=visib.m), shape=21, color="white", size=3, stroke=1) +
  geom_text( data=metar.T12.df, aes(x=lon.deg, y=lat.deg, label=radar),            color="white", hjust=0.5, vjust=2.0) +
  coord_fixed() +
  #scale_color_brewer(palette = "PRGn") +
  ggtitle("2019/02/07 @ 14:00 UTC, WRF 3-km VIS (K-G)") +
  theme(legend.position      = "bottom") +
  theme(legend.key.width     = unit(2, "cm"))

# ggplot of model Time #13 1 visib, with midwest map and color filled circles from METAR
ggplot(data = model.T13.df, aes(x=long, y=lat)) +  
  geom_point(aes(color = visib.m)) +
  geom_point(data=icicle.domain, aes(x=long, y=lat), fill=NA, color="grey50", size=0.25) +
  xlim(-98, -85) +
  ylim( 36,  50) +
  geom_point(data=metar.T13.df, aes(x=lon.deg, y=lat.deg, fill=visib.m), shape=21, color="white", size=3, stroke=1) +
  geom_text( data=metar.T13.df, aes(x=lon.deg, y=lat.deg, label=radar),            color="white", hjust=0.5, vjust=2.0) +
  coord_fixed() +
  #scale_color_brewer(palette = "PRGn") +
  ggtitle("2019/02/07 @ 15:00 UTC, WRF 3-km VIS (K-G)") +
  theme(legend.position      = "bottom") +
  theme(legend.key.width     = unit(2, "cm"))

# ggplot of model Time #14 1 visib, with midwest map and color filled circles from METAR
ggplot(data = model.T14.df, aes(x=long, y=lat)) +  
  geom_point(aes(color = visib.m)) +
  geom_point(data=icicle.domain, aes(x=long, y=lat), fill=NA, color="grey50", size=0.25) +
  xlim(-98, -85) +
  ylim( 36,  50) +
  geom_point(data=metar.T14.df, aes(x=lon.deg, y=lat.deg, fill=visib.m), shape=21, color="white", size=3, stroke=1) +
  geom_text( data=metar.T14.df, aes(x=lon.deg, y=lat.deg, label=radar),            color="white", hjust=0.5, vjust=2.0) +
  coord_fixed() +
  #scale_color_brewer(palette = "PRGn") +
  ggtitle("2019/02/07 @ 16:00 UTC, WRF 3-km VIS (K-G)") +
  theme(legend.position      = "bottom") +
  theme(legend.key.width     = unit(2, "cm"))

# ggplot of model Time #15 1 visib, with midwest map and color filled circles from METAR
ggplot(data = model.T15.df, aes(x=long, y=lat)) +  
  geom_point(aes(color = visib.m)) +
  geom_point(data=icicle.domain, aes(x=long, y=lat), fill=NA, color="grey50", size=0.25) +
  xlim(-98, -85) +
  ylim( 36,  50) +
  geom_point(data=metar.T15.df, aes(x=lon.deg, y=lat.deg, fill=visib.m), shape=21, color="white", size=3, stroke=1) +
  geom_text( data=metar.T15.df, aes(x=lon.deg, y=lat.deg, label=radar),            color="white", hjust=0.5, vjust=2.0) +
  coord_fixed() +
  #scale_color_brewer(palette = "PRGn") +
  ggtitle("2019/02/07 @ 17:00 UTC, WRF 3-km VIS (K-G)") +
  theme(legend.position      = "bottom") +
  theme(legend.key.width     = unit(2, "cm"))

# ggplot of model Time #16 1 visib, with midwest map and color filled circles from METAR
ggplot(data = model.T16.df, aes(x=long, y=lat)) +  
  geom_point(aes(color = visib.m)) +
  geom_point(data=icicle.domain, aes(x=long, y=lat), fill=NA, color="grey50", size=0.25) +
  xlim(-98, -85) +
  ylim( 36,  50) +
  geom_point(data=metar.T16.df, aes(x=lon.deg, y=lat.deg, fill=visib.m), shape=21, color="white", size=3, stroke=1) +
  geom_text( data=metar.T16.df, aes(x=lon.deg, y=lat.deg, label=radar),            color="white", hjust=0.5, vjust=2.0) +
  coord_fixed() +
  #scale_color_brewer(palette = "PRGn") +
  ggtitle("2019/02/07 @ 18:00 UTC, WRF 3-km VIS (K-G)") +
  theme(legend.position      = "bottom") +
  theme(legend.key.width     = unit(2, "cm"))

# ggplot of model Time #17 1 visib, with midwest map and color filled circles from METAR
ggplot(data = model.T17.df, aes(x=long, y=lat)) +  
  geom_point(aes(color = visib.m)) +
  geom_point(data=icicle.domain, aes(x=long, y=lat), fill=NA, color="grey50", size=0.25) +
  xlim(-98, -85) +
  ylim( 36,  50) +
  geom_point(data=metar.T17.df, aes(x=lon.deg, y=lat.deg, fill=visib.m), shape=21, color="white", size=3, stroke=1) +
  geom_text( data=metar.T17.df, aes(x=lon.deg, y=lat.deg, label=radar),            color="white", hjust=0.5, vjust=2.0) +
  coord_fixed() +
  #scale_color_brewer(palette = "PRGn") +
  ggtitle("2019/02/07 @ 19:00 UTC, WRF 3-km VIS (K-G)") +
  theme(legend.position      = "bottom") +
  theme(legend.key.width     = unit(2, "cm"))

# ggplot of model Time #18 1 visib, with midwest map and color filled circles from METAR
ggplot(data = model.T18.df, aes(x=long, y=lat)) +  
  geom_point(aes(color = visib.m)) +
  geom_point(data=icicle.domain, aes(x=long, y=lat), fill=NA, color="grey50", size=0.25) +
  xlim(-98, -85) +
  ylim( 36,  50) +
  geom_point(data=metar.T18.df, aes(x=lon.deg, y=lat.deg, fill=visib.m), shape=21, color="white", size=3, stroke=1) +
  geom_text( data=metar.T18.df, aes(x=lon.deg, y=lat.deg, label=radar),            color="white", hjust=0.5, vjust=2.0) +
  coord_fixed() +
  #scale_color_brewer(palette = "PRGn") +
  ggtitle("2019/02/07 @ 20:00 UTC, WRF 3-km VIS (K-G)") +
  theme(legend.position      = "bottom") +
  theme(legend.key.width     = unit(2, "cm"))

# ggplot of model Time #19 1 visib, with midwest map and color filled circles from METAR
ggplot(data = model.T19.df, aes(x=long, y=lat)) +  
  geom_point(aes(color = visib.m)) +
  geom_point(data=icicle.domain, aes(x=long, y=lat), fill=NA, color="grey50", size=0.25) +
  xlim(-98, -85) +
  ylim( 36,  50) +
  geom_point(data=metar.T19.df, aes(x=lon.deg, y=lat.deg, fill=visib.m), shape=21, color="white", size=3, stroke=1) +
  geom_text( data=metar.T19.df, aes(x=lon.deg, y=lat.deg, label=radar),            color="white", hjust=0.5, vjust=2.0) +
  coord_fixed() +
  #scale_color_brewer(palette = "PRGn") +
  ggtitle("2019/02/07 @ 21:00 UTC, WRF 3-km VIS (K-G)") +
  theme(legend.position      = "bottom") +
  theme(legend.key.width     = unit(2, "cm"))

# ggplot of model Time #20 1 visib, with midwest map and color filled circles from METAR
ggplot(data = model.T20.df, aes(x=long, y=lat)) +  
  geom_point(aes(color = visib.m)) +
  geom_point(data=icicle.domain, aes(x=long, y=lat), fill=NA, color="grey50", size=0.25) +
  xlim(-98, -85) +
  ylim( 36,  50) +
  geom_point(data=metar.T20.df, aes(x=lon.deg, y=lat.deg, fill=visib.m), shape=21, color="white", size=3, stroke=1) +
  geom_text( data=metar.T20.df, aes(x=lon.deg, y=lat.deg, label=radar),            color="white", hjust=0.5, vjust=2.0) +
  coord_fixed() +
  #scale_color_brewer(palette = "PRGn") +
  ggtitle("2019/02/07 @ 22:00 UTC, WRF 3-km VIS (K-G)") +
  theme(legend.position      = "bottom") +
  theme(legend.key.width     = unit(2, "cm"))


#par(mfrow=c(1,1))
#par(mar=c(4, 4, 3, 3))
#plot <- ggplot(data=data.metar.visib, aes(x=HH.MM.SS, y=visib.km*1000)) + geom_point() + geom_line() + ylim(10, 10000)
#plot + scale_y_continuous(trans='log10')

#plot(visib.metar.km*1000, type="o", pch=21, col="red", log="y", ylim=c(10, 100000), xaxt="n", xlab="Time [HH:MM]", ylab="Visibility [m]", main=paste(case.yyyymmdd, " 03:00 to 22:00 UTC @ KORD", sep=""))
#x.ticks <- obs.time.hhmm
#axis(side=1, at=x.ticks)
#axis(side=2)
#lines((data.metar.visib$visib.km[1:45]*1000), box=FALSE, type="o", pch=21, col="blue", log="y", ylim=c(10, 100000))
#legend(50, 100000, legend=c("VIS-Model", "VIS-METAR"), col=c("red", "blue"), lty=1)
#grid()
