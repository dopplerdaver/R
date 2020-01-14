#--------------------------------------------
#
# Filename: SVR.compute.VIS.SVR.dimXYZ.R
#
# Purposes: 1a. Compute SVR at every point and every height in NWP model output
#           1b. Compare SVR to RVR at many locations for case study, w/ statistical measures
#           2.  Compare model (sfc) to metar vis at many locations for case study, w/ statistical measures
#           3.  Glare/glint feasability
#
# Usage:    no initial arguements, code asks six questions of user, need to manually edit input file paths in code
#
# Inputs:   NWP model VIS and height [m, MSL] from nc file, metar vis from ? file
#
# Process:  1.  load libraries
#           2.  define constants
#           3.  solicit user input
#           4.  convert units
#           5.  define data path/filenames
#           6.  ingest input data sources
#           7.  
#           8.  plotting/ output
#
# Created:  9.20.2019 dserke
#
# THINGS TO DO:
#           1. ensure proper handling of cardinal direction VIS-H data
#           2. create/ingest multiple wrf levels (named*.L1.*, for lowest 5-10?) and calculate SVR based on user inputs (view direction, angle, etc and height diffs between model layers)
#           3. add in the 'relative airmass' = m = sec(theta), where theta is zenith angle (positive = up)
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
library(imager)
#library(rgdal)   
library(sp)
#library(rspatial) # need to load.package
library(raster)

library(ggplot2)
library(ggmap)
library(ggforce)
library(maps)
library(mapdata)

library(ncdf4)
library(NISTunits)

require(utils)
require(rgdal)
require(raster)

#--------------------
# SERVICE ACQUISITION
#   send Dave's google API account info
#register_google(key = "AIzaSyA3PxhnKvwHlHAQEdxlmVf4fvpBYD8Z0W8")     # copied directly from Google Console via 'copy' button
register_google(key = "AIzaSyBoSkGvKYWufZCQ33YcxzVfUKA89F55iw4")     # copied directly from Google Console via 'copy' button
#--------------------

#-----------------------
# define constants
#-----------------------
# algorithm parameters
num.neighbor.pixels     <- 10         # user configures, number of 3km model pixels in a given direction

# flight rule parameters
#   interpret: 'category defined from value shown and above, to value from next category.'
#   rule order (low to high vis): LIFR > IFR > MVFR > VFR
VFR.vis.nm              <- 5            # min vis (RVR) [nm] for VFR rules, below which MIFR/IFR
MVFR.vis.nm             <- 3            # [nm]
IFR.vis.nm              <- 1            # [nm]
LIFR.vis.nm             <- 0            # [nm]
VFR.ceil.ft             <- 3000         # [ft]
MVFR.ceil.ft            <- 1000         # [ft]
IFR.ceil.ft             <-  500         # [ft]
LIFR.ceil.ft            <-    0         # [ft]

# case parameters
case.yyyymmdd           <- "20190207"   # yyyymmdd
case.hh                 <- "14"         #
case.mm                 <- "00"         #

# model parameters 
NWP.model               <- "WRF"        #
NWP.horiz.res.m         <- 3000         # [m]

# location parameters
lat.KRFD.deg             <-  42.2711    # [deg], Rockford, IL
lon.KRFD.deg             <- -89.0940    # [deg]
elev.KRFD.m              <- 226         # [m]
lat.KORD.deg             <-  41.9745219 # [deg], Chicago O'Hare, IL
lon.KORD.deg             <- -87.9065972 # [deg]
elev.KORD.m              <- 207         # [m]
lat.KAAA.deg             <-  40.10      # [deg], Lincoln, IL
lon.KAAA.deg             <- -89.20      # [deg]
elev.KAAA.m              <- 182         # [m]
lat.KMWC.deg             <-  43.07      # [deg], Milwaukee/Lawren, WI
lon.KMWC.deg             <- -88.02      # [deg]
elev.KMWC.m              <- 227         # [m]
lat.KDVN.deg             <-  41.37+0.2      # [deg], Davenport, IA
lon.KDVN.deg             <- -90.35      # [deg]
elev.KDVN.m              <- 230         # [m]
lat.KDSM.deg             <-  41.32+0.2      # [deg], Des Moines, IA
lon.KDSM.deg             <- -93.40      # [deg]
elev.KDSM.m              <- 295         # [m]
lat.KAZO.deg             <-  42.14+0.2      # [deg], Kalamazoo, MI
lon.KAZO.deg             <- -85.33      # [deg]
elev.KAZO.m              <- 272         # [m]
lat.KLWA.deg             <-  42.21+0.2      # [deg], South Haven, MI
lon.KLWA.deg             <- -86.15      # [deg]
elev.KLWA.m              <- 203         # [m]
lat.K3LF.deg             <-  39.10+0.2      # [deg], Litchfield, IL
lon.K3LF.deg             <- -89.40      # [deg]
elev.K3LF.m              <- 210         # [m]
lat.KDLH.deg             <-  46.51+0.2      # [deg], Duluth, MN
lon.KDLH.deg             <- -92.12      # [deg]
elev.KDLH.m              <- 435         # [m]

#-----------------------------
# solicit user input fields
#-----------------------------

print("TO TEST SVR vs RVR ALGO IN THIS SCRIPT FOR 14Z ON 20190207 AT KORD, ENTER 0, -5, 41.97, -87.91, 730, 207")

# solicit desired viewing azimuth direction
user.input.direction         <- readline(prompt = "Q1: Choose azimuth (cardinal) direction viewer is looking [float value, either 0 for N, 90 for E, 180 for S, or 270 for W]: ")
print(paste("User choice Q1: ", user.input.direction, " deg", sep=""))

# solicit desired viewing zenith angle [deg] (neg for downward, pos for upward)
user.input.zenith.deg        <- readline(prompt = "Q2: Choose zenith angle from viewer's eye [float value from -90 to 90, deg]: ")
print(paste("User choice Q2: ", user.input.zenith.deg, " deg", sep=""))

# solicit desired viewing lat
user.input.lat.deg           <- readline(prompt = "Q3: Choose latitude of viewer's eye [float value from -90 to 90, deg]: ")
print(paste("User choice Q3: ", user.input.lat.deg, " deg", sep=""))

# solicit desired viewing lon
user.input.lon.deg           <- readline(prompt = "Q4: Choose longitude of viewer's eye [float value from -130 to -70, deg]: ")
print(paste("User choice Q4: ", user.input.lon.deg, " deg", sep=""))

# solicit desired viewing height [m, MSL]
user.input.viewer.alt.m.msl  <- readline(prompt = "Q5: Choose height of viewer's eye [float value from 0 to 15000, m MSL]: ")
print(paste("User choice Q5: ", user.input.viewer.alt.m.msl, " m", sep=""))

# solicit desired target height [m, MSL]
user.input.target.alt.m.msl  <- readline(prompt = "Q6: Choose height of viewing target [float value from 0 to 15000, m MSL]: ")
print(paste("User choice Q6: ", user.input.target.alt.m.msl, " m", sep=""))

#-----------------------------
# convert units
#-----------------------------
mPERkm                       <- 1000
VFR.vis.m                    <- NISTmileNauticalTOmeter(VFR.vis.nm)
#user.input.theta.rad         <- NISTdegTOradian(user.input.theta.deg)

#-----------------------------
# define input data file path and name
#-----------------------------
#   define WRF-3km vis/lat/lon, lowest-level-only nc file
file.path.wrf.visib          <- file.path("/d1/serke/projects/SLANT_VIS_RANGE_DETECT_FAA/data/wrf_visib/20190207/1400Z.LNUM/")
file.name.wrf.L1.visib       <- "20190207_1400.L1.nc"
file.name.wrf.L2.visib       <- "20190207_1400.L2.nc"
file.name.wrf.L3.visib       <- "20190207_1400.L3.nc"
file.name.wrf.L4.visib       <- "20190207_1400.L4.nc"
file.name.wrf.L5.visib       <- "20190207_1400.L5.nc"
file.name.wrf.L6.visib       <- "20190207_1400.L6.nc"
file.name.wrf.L7.visib       <- "20190207_1400.L7.nc"
file.name.wrf.L8.visib       <- "20190207_1400.L8.nc"
file.name.wrf.L9.visib       <- "20190207_1400.L9.nc"
file.name.wrf.L10.visib      <- "20190207_1400.L10.nc"

# define NEXRAD site location csv file path
nexrad.site.dataframetxt.dir <- file.path("/d1/serke/projects/case_studies/SNOWIE/data/RadIA_data/nexrad_site_data/")

#-----------------------------
# ingest all input data sources
#-----------------------------
#   load the NEXRAD site location text file
NEXRAD.site.df                <- read.csv(paste(nexrad.site.dataframetxt.dir, "nexrad_site.csv", sep = ""), header = FALSE, sep = ",", dec = ".", stringsAsFactors=FALSE)
colnames(NEXRAD.site.df)      <- c("NCDCID", "ICAO", "WBAN", "radname", "COUNTRY", "STATE", "COUNTY", "lat", "lon", "elev", "GMTdiff", "STN_TYPE")
head(NEXRAD.site.df)

# use model level height data from one sample time, a kluge for now but not unreasonable
model.lvl.hgts.atKORD.m       <- c(246.76785, 291.20749, 337.75629, 386.51926, 437.59088, 491.11136, 547.19183, 605.91620, 667.39801, 731.76642,  799.16156, 869.70911, 943.59399,  1021.12646,  1102.86487, 1189.14087, 1280.17859,  1375.98535, 1476.66431, 1582.34326,   1693.18872,   1809.44556,   1931.51147,   2059.72827,   2194.46558,   2336.10156,   2485.01196,   2641.55200,   2806.08008,   2978.95239,   3160.50781,   3351.16016,   3551.42725,   3761.85791,   3982.87451,   4214.54248,   4457.07178,   4710.66113,   4975.68408,   5252.60205,   5542.25391,   5845.36621,   6162.43115,   6493.38086,   6823.15625,   7151.84570,   7479.68066,   7806.93213,   8133.75732,   8460.07129,   8785.40820,   9109.31543,   9431.75488,   9754.94531,  10080.73730,  10409.13965,  10739.74805,  11089.48047,  11460.24609,  11853.43555,  12270.32617,  12713.17285,  13185.04980,  13688.41211,  14224.79102,  14796.59766,  15404.41406,  16048.11035,  16731.59766,  17462.07031,  18051.42188)
model.levl.hgts.abvSFC.m      <- model.lvl.hgts.atKORD.m - elev.KORD.m

# find the closest model height level to the viewer's height [MSL]
ind.model.lvl.closest.to.viewer<- which.min(abs(as.numeric(user.input.viewer.alt.m.msl) - model.lvl.hgts.atKORD.m))

# find the closest model height level to the target's height [MSL]
ind.model.lvl.closest.to.target<- which.min(abs(as.numeric(user.input.target.alt.m.msl) - model.lvl.hgts.atKORD.m))

num.layers.for.SVR.calc        <- abs(ind.model.lvl.closest.to.viewer - ind.model.lvl.closest.to.target)
print(paste("num.layers.for.SVR.calc = ", num.layers.for.SVR.calc, sep=""))

#   load the first 10 lowest model heights (for now)
#   ingest WRF-3km vis/lat/lon nc data, case height index #1 (MODEL.LEVEL[1] = lowest level), case time index #45 (time.model.hhmm[45] = 14:00 UTC) file
data.wrf.visib                <- nc_open(paste(file.path.wrf.visib, file.name.wrf.L1.visib, sep=""), write = FALSE, verbose = FALSE)
print(paste("The file has", data.wrf.visib$nvars, "variables"))
data.wrf.visib.var.num        <- seq(1, data.wrf.visib$nvars, by=1)
for (i in 1:length(data.wrf.visib.var.num)) {
  data.wrf.visib.nam        <- paste("v", data.wrf.visib.var.num[i], sep = "")
  assign(data.wrf.visib.nam, data.wrf.visib$var[[data.wrf.visib.var.num[i]]])
}  # end of for (i in 1:length())
if (data.wrf.visib$nvars == 4) {
  VISIB.wrf.L1.raw             <- ncvar_get( data.wrf.visib, v1 )
  XLAT.wrf.L1.raw              <- ncvar_get( data.wrf.visib, v2 )
  XLON.wrf.L1.raw              <- ncvar_get( data.wrf.visib, v3 )
  XTIME.wrf.L1.raw             <- ncvar_get( data.wrf.visib, v4 )
} else if (data.wrf.visib$nvars == 1) {
  print("Error")
}  # end of if (data.wrf.visib ...)
nc_close(data.wrf.visib)

##   ingest WRF-3km vis/lat/lon nc data, case height index #2 (MODEL.LEVEL[2] =  second lowest level), case time index #45 (time.model.hhmm[45] = 14:00 UTC) file
data.wrf.visib                <- nc_open(paste(file.path.wrf.visib, file.name.wrf.L2.visib, sep=""), write = FALSE, verbose = FALSE)
print(paste("The file has", data.wrf.visib$nvars, "variables"))
data.wrf.visib.var.num        <- seq(1, data.wrf.visib$nvars, by=1)
for (i in 1:length(data.wrf.visib.var.num)) {
  data.wrf.visib.nam        <- paste("v", data.wrf.visib.var.num[i], sep = "")
  assign(data.wrf.visib.nam, data.wrf.visib$var[[data.wrf.visib.var.num[i]]])
}  # end of for (i in 1:length())
if (data.wrf.visib$nvars == 4) {
  VISIB.wrf.L2.raw             <- ncvar_get( data.wrf.visib, v1 )
  XLAT.wrf.L2.raw              <- ncvar_get( data.wrf.visib, v2 )
  XLON.wrf.L2.raw              <- ncvar_get( data.wrf.visib, v3 )
  XTIME.wrf.L2.raw             <- ncvar_get( data.wrf.visib, v4 )
} else if (data.wrf.visib$nvars == 1) {
  print("Error")
}  # end of if (data.wrf.visib ...)
nc_close(data.wrf.visib)

##   ingest WRF-3km vis/lat/lon nc data, case height index #3 (MODEL.LEVEL[3] =  second lowest level), case time index #45 (time.model.hhmm[45] = 14:00 UTC) file
data.wrf.visib                <- nc_open(paste(file.path.wrf.visib, file.name.wrf.L3.visib, sep=""), write = FALSE, verbose = FALSE)
print(paste("The file has", data.wrf.visib$nvars, "variables"))
data.wrf.visib.var.num        <- seq(1, data.wrf.visib$nvars, by=1)
for (i in 1:length(data.wrf.visib.var.num)) {
  data.wrf.visib.nam        <- paste("v", data.wrf.visib.var.num[i], sep = "")
  assign(data.wrf.visib.nam, data.wrf.visib$var[[data.wrf.visib.var.num[i]]])
}  # end of for (i in 1:length())
if (data.wrf.visib$nvars == 4) {
  VISIB.wrf.L3.raw             <- ncvar_get( data.wrf.visib, v1 )
  XLAT.wrf.L3.raw              <- ncvar_get( data.wrf.visib, v2 )
  XLON.wrf.L3.raw              <- ncvar_get( data.wrf.visib, v3 )
  XTIME.wrf.L3.raw             <- ncvar_get( data.wrf.visib, v4 )
} else if (data.wrf.visib$nvars == 1) {
  print("Error")
}  # end of if (data.wrf.visib ...)
nc_close(data.wrf.visib)

##   ingest WRF-3km vis/lat/lon nc data, case height index #4 (MODEL.LEVEL[4] =  second lowest level), case time index #45 (time.model.hhmm[45] = 14:00 UTC) file
data.wrf.visib                <- nc_open(paste(file.path.wrf.visib, file.name.wrf.L4.visib, sep=""), write = FALSE, verbose = FALSE)
print(paste("The file has", data.wrf.visib$nvars, "variables"))
data.wrf.visib.var.num        <- seq(1, data.wrf.visib$nvars, by=1)
for (i in 1:length(data.wrf.visib.var.num)) {
  data.wrf.visib.nam        <- paste("v", data.wrf.visib.var.num[i], sep = "")
  assign(data.wrf.visib.nam, data.wrf.visib$var[[data.wrf.visib.var.num[i]]])
}  # end of for (i in 1:length())
if (data.wrf.visib$nvars == 4) {
  VISIB.wrf.L4.raw             <- ncvar_get( data.wrf.visib, v1 )
  XLAT.wrf.L4.raw              <- ncvar_get( data.wrf.visib, v2 )
  XLON.wrf.L4.raw              <- ncvar_get( data.wrf.visib, v3 )
  XTIME.wrf.L4.raw             <- ncvar_get( data.wrf.visib, v4 )
} else if (data.wrf.visib$nvars == 1) {
  print("Error")
}  # end of if (data.wrf.visib ...)
nc_close(data.wrf.visib)

##   ingest WRF-3km vis/lat/lon nc data, case height index #5 (MODEL.LEVEL[5] =  second lowest level), case time index #45 (time.model.hhmm[45] = 14:00 UTC) file
data.wrf.visib                <- nc_open(paste(file.path.wrf.visib, file.name.wrf.L5.visib, sep=""), write = FALSE, verbose = FALSE)
print(paste("The file has", data.wrf.visib$nvars, "variables"))
data.wrf.visib.var.num        <- seq(1, data.wrf.visib$nvars, by=1)
for (i in 1:length(data.wrf.visib.var.num)) {
  data.wrf.visib.nam        <- paste("v", data.wrf.visib.var.num[i], sep = "")
  assign(data.wrf.visib.nam, data.wrf.visib$var[[data.wrf.visib.var.num[i]]])
}  # end of for (i in 1:length())
if (data.wrf.visib$nvars == 4) {
  VISIB.wrf.L5.raw             <- ncvar_get( data.wrf.visib, v1 )
  XLAT.wrf.L5.raw              <- ncvar_get( data.wrf.visib, v2 )
  XLON.wrf.L5.raw              <- ncvar_get( data.wrf.visib, v3 )
  XTIME.wrf.L5.raw             <- ncvar_get( data.wrf.visib, v4 )
} else if (data.wrf.visib$nvars == 1) {
  print("Error")
}  # end of if (data.wrf.visib ...)
nc_close(data.wrf.visib)

##   ingest WRF-3km vis/lat/lon nc data, case height index #6 (MODEL.LEVEL[6] =  second lowest level), case time index #45 (time.model.hhmm[45] = 14:00 UTC) file
data.wrf.visib                <- nc_open(paste(file.path.wrf.visib, file.name.wrf.L6.visib, sep=""), write = FALSE, verbose = FALSE)
print(paste("The file has", data.wrf.visib$nvars, "variables"))
data.wrf.visib.var.num        <- seq(1, data.wrf.visib$nvars, by=1)
for (i in 1:length(data.wrf.visib.var.num)) {
  data.wrf.visib.nam        <- paste("v", data.wrf.visib.var.num[i], sep = "")
  assign(data.wrf.visib.nam, data.wrf.visib$var[[data.wrf.visib.var.num[i]]])
}  # end of for (i in 1:length())
if (data.wrf.visib$nvars == 4) {
  VISIB.wrf.L6.raw             <- ncvar_get( data.wrf.visib, v1 )
  XLAT.wrf.L6.raw              <- ncvar_get( data.wrf.visib, v2 )
  XLON.wrf.L6.raw              <- ncvar_get( data.wrf.visib, v3 )
  XTIME.wrf.L6.raw             <- ncvar_get( data.wrf.visib, v4 )
} else if (data.wrf.visib$nvars == 1) {
  print("Error")
}  # end of if (data.wrf.visib ...)
nc_close(data.wrf.visib)

##   ingest WRF-3km vis/lat/lon nc data, case height index #7 (MODEL.LEVEL[7] =  second lowest level), case time index #45 (time.model.hhmm[45] = 14:00 UTC) file
data.wrf.visib                <- nc_open(paste(file.path.wrf.visib, file.name.wrf.L7.visib, sep=""), write = FALSE, verbose = FALSE)
print(paste("The file has", data.wrf.visib$nvars, "variables"))
data.wrf.visib.var.num        <- seq(1, data.wrf.visib$nvars, by=1)
for (i in 1:length(data.wrf.visib.var.num)) {
  data.wrf.visib.nam        <- paste("v", data.wrf.visib.var.num[i], sep = "")
  assign(data.wrf.visib.nam, data.wrf.visib$var[[data.wrf.visib.var.num[i]]])
}  # end of for (i in 1:length())
if (data.wrf.visib$nvars == 4) {
  VISIB.wrf.L7.raw             <- ncvar_get( data.wrf.visib, v1 )
  XLAT.wrf.L7.raw              <- ncvar_get( data.wrf.visib, v2 )
  XLON.wrf.L7.raw              <- ncvar_get( data.wrf.visib, v3 )
  XTIME.wrf.L7.raw             <- ncvar_get( data.wrf.visib, v4 )
} else if (data.wrf.visib$nvars == 1) {
  print("Error")
}  # end of if (data.wrf.visib ...)
nc_close(data.wrf.visib)

##   ingest WRF-3km vis/lat/lon nc data, case height index #8 (MODEL.LEVEL[8] =  second lowest level), case time index #45 (time.model.hhmm[45] = 14:00 UTC) file
data.wrf.visib                <- nc_open(paste(file.path.wrf.visib, file.name.wrf.L8.visib, sep=""), write = FALSE, verbose = FALSE)
print(paste("The file has", data.wrf.visib$nvars, "variables"))
data.wrf.visib.var.num        <- seq(1, data.wrf.visib$nvars, by=1)
for (i in 1:length(data.wrf.visib.var.num)) {
  data.wrf.visib.nam        <- paste("v", data.wrf.visib.var.num[i], sep = "")
  assign(data.wrf.visib.nam, data.wrf.visib$var[[data.wrf.visib.var.num[i]]])
}  # end of for (i in 1:length())
if (data.wrf.visib$nvars == 4) {
  VISIB.wrf.L8.raw             <- ncvar_get( data.wrf.visib, v1 )
  XLAT.wrf.L8.raw              <- ncvar_get( data.wrf.visib, v2 )
  XLON.wrf.L8.raw              <- ncvar_get( data.wrf.visib, v3 )
  XTIME.wrf.L8.raw             <- ncvar_get( data.wrf.visib, v4 )
} else if (data.wrf.visib$nvars == 1) {
  print("Error")
}  # end of if (data.wrf.visib ...)
nc_close(data.wrf.visib)

##   ingest WRF-3km vis/lat/lon nc data, case height index #9 (MODEL.LEVEL[9] =  second lowest level), case time index #45 (time.model.hhmm[45] = 14:00 UTC) file
data.wrf.visib                <- nc_open(paste(file.path.wrf.visib, file.name.wrf.L9.visib, sep=""), write = FALSE, verbose = FALSE)
print(paste("The file has", data.wrf.visib$nvars, "variables"))
data.wrf.visib.var.num        <- seq(1, data.wrf.visib$nvars, by=1)
for (i in 1:length(data.wrf.visib.var.num)) {
  data.wrf.visib.nam        <- paste("v", data.wrf.visib.var.num[i], sep = "")
  assign(data.wrf.visib.nam, data.wrf.visib$var[[data.wrf.visib.var.num[i]]])
}  # end of for (i in 1:length())
if (data.wrf.visib$nvars == 4) {
  VISIB.wrf.L9.raw             <- ncvar_get( data.wrf.visib, v1 )
  XLAT.wrf.L9.raw              <- ncvar_get( data.wrf.visib, v2 )
  XLON.wrf.L9.raw              <- ncvar_get( data.wrf.visib, v3 )
  XTIME.wrf.L9.raw             <- ncvar_get( data.wrf.visib, v4 )
} else if (data.wrf.visib$nvars == 1) {
  print("Error")
}  # end of if (data.wrf.visib ...)
nc_close(data.wrf.visib)

##   ingest WRF-3km vis/lat/lon nc data, case height index #10 (MODEL.LEVEL[10] =  second lowest level), case time index #45 (time.model.hhmm[45] = 14:00 UTC) file
data.wrf.visib                <- nc_open(paste(file.path.wrf.visib, file.name.wrf.L10.visib, sep=""), write = FALSE, verbose = FALSE)
print(paste("The file has", data.wrf.visib$nvars, "variables"))
data.wrf.visib.var.num        <- seq(1, data.wrf.visib$nvars, by=1)
for (i in 1:length(data.wrf.visib.var.num)) {
  data.wrf.visib.nam        <- paste("v", data.wrf.visib.var.num[i], sep = "")
  assign(data.wrf.visib.nam, data.wrf.visib$var[[data.wrf.visib.var.num[i]]])
}  # end of for (i in 1:length())
if (data.wrf.visib$nvars == 4) {
  VISIB.wrf.L10.raw             <- ncvar_get( data.wrf.visib, v1 )
  XLAT.wrf.L10.raw              <- ncvar_get( data.wrf.visib, v2 )
  XLON.wrf.L10.raw              <- ncvar_get( data.wrf.visib, v3 )
  XTIME.wrf.L10.raw             <- ncvar_get( data.wrf.visib, v4 )
} else if (data.wrf.visib$nvars == 1) {
  print("Error")
}  # end of if (data.wrf.visib ...)
nc_close(data.wrf.visib)

#-----------------------------
# create VISIB 2-D matrix that recomputes vis in multiple directions accounting for user configured number of neighbor values
#-----------------------------
# demonstrate capability by computing for only N, E, S, and W of pixel[ii, jj]
# NOTE: VISIB.raw[1, 1] is SW corner

# initialize new VISIB 2-D matrices
VISIB.wrf.L1.toN.neighbors          <- zeros(dim(VISIB.wrf.L1.raw)[1], dim(VISIB.wrf.L1.raw)[2])
VISIB.wrf.L2.toN.neighbors          <- zeros(dim(VISIB.wrf.L2.raw)[1], dim(VISIB.wrf.L2.raw)[2])
VISIB.wrf.L3.toN.neighbors          <- zeros(dim(VISIB.wrf.L3.raw)[1], dim(VISIB.wrf.L3.raw)[2])
VISIB.wrf.L4.toN.neighbors          <- zeros(dim(VISIB.wrf.L4.raw)[1], dim(VISIB.wrf.L4.raw)[2])
VISIB.wrf.L5.toN.neighbors          <- zeros(dim(VISIB.wrf.L5.raw)[1], dim(VISIB.wrf.L5.raw)[2])
VISIB.wrf.L6.toN.neighbors          <- zeros(dim(VISIB.wrf.L6.raw)[1], dim(VISIB.wrf.L6.raw)[2])
VISIB.wrf.L7.toN.neighbors          <- zeros(dim(VISIB.wrf.L7.raw)[1], dim(VISIB.wrf.L7.raw)[2])
VISIB.wrf.L8.toN.neighbors          <- zeros(dim(VISIB.wrf.L8.raw)[1], dim(VISIB.wrf.L8.raw)[2])
VISIB.wrf.L9.toN.neighbors          <- zeros(dim(VISIB.wrf.L9.raw)[1], dim(VISIB.wrf.L9.raw)[2])
VISIB.wrf.L10.toN.neighbors         <- zeros(dim(VISIB.wrf.L10.raw)[1], dim(VISIB.wrf.L10.raw)[2])
#VISIB.wrf.L1.toE.neighbors          <- zeros(dim(VISIB.wrf.L1.raw)[1], dim(VISIB.wrf.L1.raw)[2])
#VISIB.wrf.L1.toS.neighbors          <- zeros(dim(VISIB.wrf.L1.raw)[1], dim(VISIB.wrf.L1.raw)[2])
#VISIB.wrf.L1.toW.neighbors          <- zeros(dim(VISIB.wrf.L1.raw)[1], dim(VISIB.wrf.L1.raw)[2])

# compute facing N of pixel[ii, jj]
# loop over all lon (ii, interior first) and lat (jj, exterior second)
for (jj in 1:dim(VISIB.wrf.L1.raw)[2]) {
  
  for (ii in 1:dim(VISIB.wrf.L1.raw)[1]) {
    
    if (jj < dim(VISIB.wrf.L1.raw)[2] - num.neighbor.pixels) {
      
      for (kk in 1:num.neighbor.pixels) {
        
        # if vis decreasing in this direction then set vis to min vis that is closest
        if (VISIB.wrf.L1.raw[ii, jj] < VISIB.wrf.L1.raw[ii, jj+kk]) {
          VISIB.wrf.L1.toN.neighbors[ii, jj] <- VISIB.wrf.L1.raw[ii, jj]
          kk                          <- num.neighbor.pixels # to terminate kk for loop
          
        # else if gt horiz model resolution then do nothing and move to next further away pixel
        } else if (VISIB.wrf.L1.raw[ii, jj] >= VISIB.wrf.L1.raw[ii, jj+kk] & VISIB.wrf.L1.raw[ii, jj+kk] > NWP.horiz.res.m) {
          # move on to next neighboring pixel vis value
          kk                          <- kk
        
        # else further out pixel vis lt horiz model resolution then increment vis by horiz model rez plus neighbor pixel res until at max num.neighbor.pixels
        } else {
          VISIB.wrf.L1.toN.neighbors[ii, jj] <- VISIB.wrf.L1.raw[ii, jj+kk] + (NWP.horiz.res.m * kk) # in meters
          kk                          <- num.neighbor.pixels # to terminate kk for loop
        }  # end of if (VISIB.wrf.raw[ii, jj] ...)
      
      }    # end of for (kk in 1:num.neighbor.pixels) {}
      
    }      # end of if (jj < dim(VISIB.wrf.raw)[1]) {}
    
  }        # end of for (ii in ...)
  
}          # end of for (jj in ...)

## compute facing E of pixel[ii, jj]
## loop over all lon (ii, interior first) and lat (jj, exterior second)
#for (jj in 1:dim(VISIB.wrf.raw)[2]) {
#  
#  for (ii in 1:dim(VISIB.wrf.raw)[1]) {
#
#    if (ii < dim(VISIB.wrf.raw)[2] - num.neighbor.pixels) {
#  
#      for (kk in 1:num.neighbor.pixels) {
#    
#        # if vis decreasing in this direction then set vis to min vis that is closest
#        if (VISIB.wrf.raw[ii, jj] < VISIB.wrf.raw[ii+kk, jj]) {
#          VISIB.wrf.toN.neighbors[ii, jj] <- VISIB.wrf.raw[ii, jj]
#          kk                          <- num.neighbor.pixels # to terminate kk for loop
#      
#        # else if gt horiz model resolution then do nothing and move to next further away pixel
#        } else if (VISIB.wrf.raw[ii, jj] >= VISIB.wrf.raw[ii+kk, jj] & VISIB.wrf.raw[ii+kk, jj] > NWP.horiz.res.m) {
#          # move on to next neighboring pixel vis value
#          kk                          <- kk
#      
#        # else further out pixel vis lt horiz model resolution then increment vis by horiz model rez plus neighbor pixel res until at max num.neighbor.pixels
#        } else {
#          VISIB.wrf.toN.neighbors[ii, jj] <- VISIB.wrf.raw[ii+kk, jj] + (NWP.horiz.res.m * kk) # in meters
#          kk                          <- num.neighbor.pixels # to terminate kk for loop
#        }  # end of if (VISIB.wrf.raw[ii, jj] ...)
#    
#      }    # end of for (kk in 1:num.neighbor.pixels) {}
#  
#    }      # end of if (ii < dim(VISIB.wrf.raw)[1]) {}
#
#  }  # end of for (ii in ...)
#
#}    # end of for (jj in ...)

# NEED TO COPY AND EDIT LOOP CHUNKS ABOVE FOR E, S AND W CARDINAL DIRECTIONS

#-----------------------------
# derive SVR in a given azimuthal direction from columnar vis-H values
#-----------------------------
# NOTE: SVR must be calculated at least one model layer in depth (depth = opposite side)
# TRIG: Cos(theta.rad) = (adjacent side) / (hypotenuse side) = VISIB.wrf.to{dir}.m / SVR.to{dir}.m
#       - user provides theta.deg near top, converted to theta.rad
#       SOLVE FOR SVR.to{dir}.m = VISIB.wrf.to{dir}.m / Cos(theta.rad)
# SVR(downward).to{dir} = layer.depth.m
#SVR.toN.m                          <- VISIB.wrf.toN.neighbors / cos(user.input.theta.rad)
#SVR.toE.m                          <- VISIB.wrf.toE.neighbors / cos(user.input.theta.rad)
#SVR.toS.m                          <- VISIB.wrf.toS.neighbors / cos(user.input.theta.rad)
#SVR.toW.m                          <- VISIB.wrf.toW.neighbors / cos(user.input.theta.rad)

# NEED TO DECIDE HOW TO DEAL WITH CALCULATING SVR SLANT DOWN THROUGH MULTIPLE MODEL LAYERS HERE

#-----------------------------
# restrict VISIB.raw to less than VFR RVR minimum [m], mainly to restrict the dynamic range of subsequent plots
#-----------------------------
ind.L1.raw.gtVFR                             <- which(VISIB.wrf.L1.raw > VFR.vis.m)
VISIB.wrf.L1.raw[ind.L1.raw.gtVFR]           <- VFR.vis.m
ind.L1.toN.gtVFR                             <- which(VISIB.wrf.L1.toN.neighbors > VFR.vis.m)
VISIB.wrf.L1.toN.neighbors[ind.L1.toN.gtVFR] <- VFR.vis.m
#ind.L1.toE.gtVFR                             <- which(VISIB.wrf.L1.toE.neighbors > VFR.vis.m)
#VISIB.wrf.L1.toE.neighbors[ind.L1.toE.gtVFR] <- VFR.vis.m
#ind.L1.toS.gtVFR                             <- which(VISIB.wrf.L1.toS.neighbors > VFR.vis.m)
#VISIB.wrf.L1.toS.neighbors[ind.L1.toS.gtVFR] <- VFR.vis.m
#ind.L1.toW.gtVFR                             <- which(VISIB.wrf.L1.toW.neighbors > VFR.vis.m)
#VISIB.wrf.L1.toW.neighbors[ind.L1.toW.gtVFR] <- VFR.vis.m

#-----------------------------
# vis/lat/lon into data frame
#-----------------------------
model.L1.df                                  <- data.frame(as.vector(VISIB.wrf.L1.raw), as.vector(XLON.wrf.L1.raw), as.vector(XLAT.wrf.L1.raw))
colnames(model.L1.df)                        <- c("visib.m", "long", "lat")
head(model.L1.df)

#-----------------------------
# build raster layer(s) and edit associated values
#-----------------------------
# convert data matrices to raster format
VISIB.wrf.L1.raw.raster                     <- raster(VISIB.wrf.L1.raw)
VISIB.wrf.L1.toN.neighbors.raster           <- raster(VISIB.wrf.L1.toN.neighbors)
#VISIB.wrf.L1.toE.neighbors.raster           <- raster(VISIB.wrf.L1.toE.neighbors)
#VISIB.wrf.L1.toS.neighbors.raster           <- raster(VISIB.wrf.L1.toS.neighbors)
#VISIB.wrf.L1.toW.neighbors.raster           <- raster(VISIB.wrf.L1.toW.neighbors)

# print some raster values
VISIB.wrf.L1.raw.raster
projection(VISIB.wrf.L1.raw.raster)

# define x/y extent
extent(VISIB.wrf.L1.raw.raster)             <- extent(min(XLAT.wrf.L1.raw), max(XLAT.wrf.L1.raw), min(XLON.wrf.L1.raw), max(XLON.wrf.L1.raw))
extent(VISIB.wrf.L1.toN.neighbors.raster)   <- extent(min(XLAT.wrf.L1.raw), max(XLAT.wrf.L1.raw), min(XLON.wrf.L1.raw), max(XLON.wrf.L1.raw))
#extent(VISIB.wrf.L1.toE.neighbors.raster)   <- extent(min(XLAT.wrf.L1.raw), max(XLAT.wrf.L1.raw), min(XLON.wrf.L1.raw), max(XLON.wrf.L1.raw))
#extent(VISIB.wrf.L1.toS.neighbors.raster)   <- extent(min(XLAT.wrf.L1.raw), max(XLAT.wrf.L1.raw), min(XLON.wrf.L1.raw), max(XLON.wrf.L1.raw))
#extent(VISIB.wrf.L1.toW.neighbors.raster)   <- extent(min(XLAT.wrf.L1.raw), max(XLAT.wrf.L1.raw), min(XLON.wrf.L1.raw), max(XLON.wrf.L1.raw))

# set raster CRS values
crs(VISIB.wrf.L1.raw.raster)                <- "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"

# print some raster values
VISIB.wrf.L1.raw.raster

# go the ggmap route for overlay raster layer spatial data with map
VISIB.wrf.L1.raw.agg.raster                 <- aggregate(VISIB.wrf.L1.raw.raster, fact=4)
rtp.L1                                      <- rasterToPolygons(VISIB.wrf.L1.raw.agg.raster)
ICICLE.map                                  <- get_map(location = c(-90, 42.5), maptype = "terrain", source = "google", zoom = 7)
ggmap(ICICLE.map)

#ICICLE.map + geom_raster(data=rtp, aes(fill=layer), alpha=0.5) +
#             scale_fill_gradientn("RasterValues", colors=topo.colors(255))

#ggmap(ICICLE.map, base_layer=ggplot(model.df, aes(x=long, y=lat, z=visib.m))) +
#  geom_raster(aes(fill=visib.m), alpha=0.5) +
#  coord_cartesian() +
#  scale_fill_gradient(low="white", high="blue") +
#  labs(fill="", title="test", x="Longitude [deg]", y="Latitude [deg]")

#-----------------------------
# plotting
#-----------------------------
# plot vis-H derived from model
par(mfrow=c(1,1))
par(mar=c(4, 4, 3, 3))
plot(t(flip((VISIB.wrf.L1.raw.raster), 1)) / mPERkm, box=FALSE, xlab="Lon [deg]", ylab="Lat [deg]", xlim=c(-97, -85), ylim=c(36, 48), main=paste(case.yyyymmdd, " ", case.hh, ":", case.mm, " UTC Vis-H(", NWP.model, "-3km @ ~40m alt AGL) < VFR", sep=""))
states                             <- map_data("state")
icicle.domain                      <- subset(states, region %in% c("minnesota", "iowa", "missouri", "kentucky", "illinois", "wisconsin", "michigan", "indiana", "ohio"))
lines(icicle.domain, add=TRUE, type="p", pch=16, col="black", cex=0.3)
grid()
#out <- draw_circle(p1, lon.KRFD.deg, lat.KRFD.deg, 5, color = c(0,1,0), opacity = 1, filled = TRUE)
#plot(out)
text(lon.KRFD.deg, lat.KRFD.deg, "*")
text(lon.KRFD.deg, lat.KRFD.deg-0.2, "KRFD")
text(lon.KORD.deg, lat.KORD.deg, "*")
text(lon.KORD.deg, lat.KORD.deg-0.2, "KORD")
text(lon.KAAA.deg, lat.KAAA.deg, "*")
text(lon.KAAA.deg, lat.KAAA.deg-0.2, "KAAA")
text(lon.KMWC.deg, lat.KMWC.deg, "*")
text(lon.KMWC.deg, lat.KMWC.deg-0.2, "KMWC")
text(lon.KDVN.deg, lat.KDVN.deg, "*")
text(lon.KDVN.deg, lat.KDVN.deg-0.2, "KDVN")
text(lon.KDSM.deg, lat.KDSM.deg, "*")
text(lon.KDSM.deg, lat.KDSM.deg-0.2, "KDSM")
text(lon.KAZO.deg, lat.KAZO.deg, "*")
text(lon.KAZO.deg, lat.KAZO.deg-0.2, "KAZO")
text(lon.KLWA.deg, lat.KLWA.deg, "*")
text(lon.KLWA.deg, lat.KLWA.deg-0.2, "KLWA")
text(lon.K3LF.deg, lat.K3LF.deg, "*")
text(lon.K3LF.deg, lat.K3LF.deg-0.2, "K3LF")
text(lon.KDLH.deg, lat.KDLH.deg, "*")
text(lon.KDLH.deg, lat.KDLH.deg-0.2, "KDLH")

# build data arrays that make up metar.df
visib.km           <- c(         4.8,          0.4,          2.4,          2.4,          0.4,          4.8,          4.8,          2.8,         11.3)
visib.m            <- visib.km * 1000
time.HHMM          <- c(     "14:07",      "14:04",      "14:07",      "13:54",      "13:55",      "13:55",      "13:55",      "13:55",      "13:55")
radar              <- c(      "KORD",       "KMWC",       "KDVN",       "KDSM",       "KAAA",       "KDLH",       "KAZO",       "KLWA",       "K3LF")
lon.deg            <- c(lon.KORD.deg, lon.KMWC.deg, lon.KDVN.deg, lon.KDSM.deg, lon.KAAA.deg, lon.KDLH.deg, lon.KAZO.deg, lon.KLWA.deg, lon.K3LF.deg)
lat.deg            <- c(lat.KORD.deg, lat.KMWC.deg, lat.KDVN.deg, lat.KDSM.deg, lat.KAAA.deg, lat.KDLH.deg, lat.KAZO.deg, lat.KLWA.deg, lat.K3LF.deg)
# build metar df, over 9 locations x 20 hours
metar.df           <- data.frame(radar, lon.deg, lat.deg, time.HHMM, visib.m)
colnames(metar.df) <- c("radar", "lon.deg", "lat.deg", "HH:MM", "visib.m")
head(metar.df)

# ggplot of model Level 1 visib, with midwest map and color filled circles from METAR
ggplot(data = model.L1.df, aes(x=long, y=lat)) +  
  geom_point(aes(color = visib.m)) +
  #geom_tile() +
  #geom_tile(aes(fill=visib.m), alpha=0.8) +
  #geom_tile(  data=model.L1.df, aes(x=long, y=lat, fill=visib.m), alpha=0.8) +
  geom_polygon(data=icicle.domain, aes(x=long, y=lat, group=group), fill=NA, color="grey50", size=0.25) +
  xlim(-98, -83) +
  ylim( 35,  50) +
  #coord_equal()  +
  geom_point(data=metar.df, aes(x=lon.deg, y=lat.deg, fill=visib.m), shape=21, color="white", size=3, stroke=1) +
  geom_text( data=metar.df, aes(x=lon.deg, y=lat.deg, label=radar),                  color="white", hjust=0.5, vjust=2.0) +
  coord_fixed() +
  #scale_color_brewer(palette = "PRGn") +
  ggtitle("2019/02/07 @ 14:00 UTC, WRF 3-km VIS (K-G)") +
  theme(legend.position      = "bottom") +
  theme(legend.key.width     = unit(2, "cm"))

#plot(t(flip((VISIB.wrf.L1.raw.raster - VISIB.wrf.L1.toN.neighbors.raster), 1)) / mPERkm, box=FALSE, xlab="Lon [deg]", ylab="Lat [deg]", main=paste(case.yyyymmdd, " ", case.hh, ":", case.mm, " UTC Vis-H(", NWP.model, "-3km) < VFR @ SFC, diff(raw-toE)", sep=""))
#grid()
#text(lon.KRFD.deg, lat.KRFD.deg, "*")
#text(lon.KORD.deg, lat.KORD.deg, "*")


# plot obs (metar) versus modelled vis at lowest medel level (sfc)


# plot derived SVR


# plot derived SVR versus RVR


# plot mapping samples
## this map plotting works but can't get raster data to plot over.  save for future reference
#states                             <- map_data("state")
#icicle.domain                      <- subset(states, region %in% c("illinois", "ohio", "iowa", "wisconsin", "michigan", "missouri", "indiana", "kentucky", "minnesota"))
#locations                          <- data.frame(long=c(lon.KRFD.deg, lon.KORD.deg), lat=c(lat.KRFD.deg, lat.KORD.deg), names=c("KRFD", "KORD"), stringsAsFactors=FALSE)
#gg1                                <- ggplot(model.df) + geom_polygon(data=icicle.domain, aes(x=long, y=lat, group=group), fill=NA, color="black") + coord_fixed(1.3)
#gg1 + geom_point(data=NEXRAD.site.df, aes(x=lon, y=lat), color="black", size=3) +
#      geom_point(data=NEXRAD.site.df, aes(x=lon, y=lat), color="yellow", size=2) +
#      geom_point(data=locations, aes(x=long[1], y=lat[1]), color="black", size=3) +
#      geom_point(data=locations, aes(x=long[1], y=lat[1]), color="blue", size=2) +
#      xlim(-97.5, -81.5) +
#      ylim(35.5, 49.5)

