#--------------------------------------------
#
# Filename: HAIL.ACCUM.detect.R
#
# Usage:    no initial arguements, edit input file path in code
#
# Purpose:  Detection of realtime and accumulated hail depth operational NEXRAD Ze
#           This algo introduces two new features:
#           1. a dynamic diameter fall velocity Vt that depends on hailstone diameter D, and
#           2. a statistically best coefficient (epsilon) which results in min diff between observed and radar-based hail accum
#           3. Will et al. (1998) method used Z-to-hail (or Z-to-Energy) relation, instead of Z-to-liquid water relation since energy relates to damage
#
# Inputs:   1. L2 Z (dBZ) from single radar
#           2. temp profile, nearest spatially/temporally,
#           3. L3 HCA Product - 6 tilt levels, specifically hail (9) and hail with rain (10) category
#           4. L3 Storm Structure Product - overlay
#
# Process:  1.  retrieve sounding from obs or model (closest time/space) and store to archive
#           2a. retrieve storm cell location, speed and direction from NEXRAD L3 SSTRUC output at: 
#                   ftp://ftp.ncdc.noaa.gov/pub/data/swdi/database-csv/v2/ filename = structure{YYYY}.csv.gz
#           2.  retrieve HCA hail and hail mixed with rain (9&10) from NEXRAD L3 HCA output at:
#                   https://www.ncdc.noaa.gov/nexradinv/chooseday.jsp?id=kama
#                   Then convert these NCDC binary format files to mdv then to netcdf format as per Sue Dettling's two email streams during 9/2019
#           2b. compute IWC.hail (eqn 2) below 0 deg to id hail/rain mix
#           2c. calculate MESH, necessary to compute Vel.t.ofdiam
#           3.  compute hail depth for this radar volume
#           4.  sum hail accum depth across all previous times steps for each sfc pixel
#           5.  produce an accum.sfc.hail, storm speed, real-time hail size maps
#
# Created:  9.20.2019 dserke
#
# THINGS TO DO:
# 1. get process running on multiple KAMA volumes
# 2. check/test output of per vol accum and tot accum for 1.5 hr test period
# 3. 
#
#
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
library(NISTunits)
library(lattice)
library(ggplot2)
library(ggmap)
library(maps)
library(mapdata)
library(geosphere)

#--------------------
# SERVICE ACQUISITION
#   send Dave's google API account info
#register_google(key = "AIzaSyA3PxhnKvwHlHAQEdxlmVf4fvpBYD8Z0W8")     # copied directly from Google Console via 'copy' button
register_google(key = "AIzaSyBoSkGvKYWufZCQ33YcxzVfUKA89F55iw4")     # copied directly from Google Console via 'copy' button
#--------------------

#-----------------------
# define constants
#-----------------------
options(max.print=10000)
flag.first.time               <- 0                       # binary flag
flag.display.to.gui           <- 1
flag.make.plots               <- 0                       # binary flag

secPERmin                     <- 60                      # time conversion factor
m3PERcm3                      <- 1e-6                    # volume conversion factor
kmTOm                         <- 1000                    # distance conversion factor
cmPERm                        <- 100                     # distance conversion factor
mPERcm                        <- 0.01
radius.earth.km               <- 6371                    # [km]
refract.index                 <- 1.21                    # [NA]

packing.density               <- 0.64                    # [na], closest possible random packing of monodisperse spheres (Scott and Kilgour, 1969)
rho.hail                      <- 0.90                    # density [g/cm3]
epsilon                       <- 0.81                    # [na], coeff for min diff between obs and radar-based hail accum
Z.subL                        <- 40                      # [dBZ] hardcoded from Witt et al paper
Z.subU                        <- 50                      # [dBZ] hardcoded from Witt et al paper
Z.value.for.echo.top          <- 18

#-----------------------------
# case study information
#-----------------------------
# https://k1025.com/vicksburg-hailstorm-august-2019/
# https://wkzo.com/news/articles/2019/aug/15/summer-hail-hits-southwest-michigan/927864/
# https://maps.hailstrike.com/tag/vicksburg-mi/
# https://wwmt.com/news/local/clean-up-continues-after-storms-damage-trees-power-lines-and-cars-in-vicksburg
# Wind Report at 3:30 P 08/14/2019 0 mi E
# SEVERAL LARGE TREES DOWN IN VICKSBURG
NEXRAD.radar                     <- "KAMA"                  # Amarillo, TX NEXRAD
case.yyyymmdd                    <- "20120411"              # [YYYYMMDD]
case.hhmm.utc                    <- c("04:30", "04:34", "04:39", "04:43", "04:47", "04:52", "04:56", "05:00", "05:05", "05:09", "05:13", "05:18", "05:22", "05:26", "05:30", "05:35")              # [UTC]
case.num.radar.volumes           <- length(case.hhmm.utc)
theta.radar.volume.L2.deg        <- c(0.5, 0.9, 1.4, 1.9, 2.4, 3.1, 3.9, 5.1, 6.4, 8.1, 10.1, 12.5, 15.4, 19.5)
theta.radar.volume.L3.deg        <- c(0.5, 0.9, 1.5, 1.8, 2.4, 3.4)

#-----------------------------
# define data file paths and names
#-----------------------------
# define NEXRAD site location csv file path
nexrad.site.dataframetxt.dir     <- file.path("/d1/serke/projects/case_studies/SNOWIE/data/RadIA_data/nexrad_site_data/")

# define balloon sounding csv file path and name
balloon.sounding.file.dir        <- file.path("/d1/serke/projects/HAIL_ACCUM_detect_ATEC/data/sndg/20120411_KAMA/")
balloon.sounding.listfiles       <- list.files(path = balloon.sounding.file.dir , include.dirs = FALSE)
balloon.sounding.file.name       <- balloon.sounding.listfiles

# define NEXRAD L2 nc file path and name
nexrad.L2.file.dir               <- file.path(paste("/d1/serke/projects/HAIL_ACCUM_detect_ATEC/data/NEXRAD_L2/", case.yyyymmdd, "_", NEXRAD.radar, "/", sep=""))
nexrad.L2.listfiles              <- list.files(path = nexrad.L2.file.dir , include.dirs = FALSE)

# calculate number of radar data volumes available
count.volumes                    <- 0
tilt.last                        <- 0 
for (k in 2:length(nexrad.L2.listfiles) - 1) {
  
  if (substr(nexrad.L2.listfiles[k], 35, 35) == "_") {
    tilt.curr <- as.numeric(substr(nexrad.L2.listfiles[k], 32, 34))
  } else {
    tilt.curr <- as.numeric(substr(nexrad.L2.listfiles[k], 32, 35))
  }  # end of if (substr(nexrad.L2.listfiles) - 1) ...
  
  if (tilt.curr < tilt.last) {
    count.volumes <- count.volumes + 1
  }
  tilt.last <- tilt.curr
}  # end of for (k in ...)

# calculate number of radar data files in each volume
tilt.last                        <- 0 
count.volumes                    <- 0
count.files.in.volume            <- 0
count.files.in.volume.array      <- array(data = 0L, dim = c(count.volumes))
for (l in 2:length(nexrad.L2.listfiles) - 1) {
  if (substr(nexrad.L2.listfiles[l], 35, 35) == "_") {
    tilt.curr <- as.numeric(substr(nexrad.L2.listfiles[l], 32, 34))
  } else {
    tilt.curr <- as.numeric(substr(nexrad.L2.listfiles[l], 32, 35))
  }  # end of if (substr(nexrad.L2.listfiles) - 1) ...
  
  if (tilt.curr < tilt.last) {
    count.volumes                              <- count.volumes + 1
    count.files.in.volume.array[count.volumes] <- count.files.in.volume
    count.files.in.volume                      <- 0
  } # end of if (tilt.curr < ...)
  
  count.files.in.volume <- count.files.in.volume + 1
  tilt.last             <- tilt.curr
  
}   # end of for (l in 2: length(nexrad.L2.listfiles[l])...)

# define NEXRAD L3 HCA data file path and names
#   file path
nexrad.L3.HCA.0.5deg.file.dir    <- file.path(paste("/d1/serke/projects/HAIL_ACCUM_detect_ATEC/data/NEXRAD_L3/", case.yyyymmdd, "_", NEXRAD.radar, "/nc/N0H/", sep=""))
nexrad.L3.HCA.0.9deg.file.dir    <- file.path(paste("/d1/serke/projects/HAIL_ACCUM_detect_ATEC/data/NEXRAD_L3/", case.yyyymmdd, "_", NEXRAD.radar, "/nc/NAH/", sep=""))
nexrad.L3.HCA.1.5deg.file.dir    <- file.path(paste("/d1/serke/projects/HAIL_ACCUM_detect_ATEC/data/NEXRAD_L3/", case.yyyymmdd, "_", NEXRAD.radar, "/nc/N1H/", sep=""))
nexrad.L3.HCA.1.8deg.file.dir    <- file.path(paste("/d1/serke/projects/HAIL_ACCUM_detect_ATEC/data/NEXRAD_L3/", case.yyyymmdd, "_", NEXRAD.radar, "/nc/NBH/", sep=""))
nexrad.L3.HCA.2.4deg.file.dir    <- file.path(paste("/d1/serke/projects/HAIL_ACCUM_detect_ATEC/data/NEXRAD_L3/", case.yyyymmdd, "_", NEXRAD.radar, "/nc/N2H/", sep=""))
nexrad.L3.HCA.3.4deg.file.dir    <- file.path(paste("/d1/serke/projects/HAIL_ACCUM_detect_ATEC/data/NEXRAD_L3/", case.yyyymmdd, "_", NEXRAD.radar, "/nc/N3H/", sep=""))
#   file names
nexrad.L3.HCA.0.5deg.listfiles   <- list.files(path = nexrad.L3.HCA.0.5deg.file.dir , include.dirs = FALSE)
nexrad.L3.HCA.0.9deg.listfiles   <- list.files(path = nexrad.L3.HCA.0.9deg.file.dir , include.dirs = FALSE)
nexrad.L3.HCA.1.5deg.listfiles   <- list.files(path = nexrad.L3.HCA.1.5deg.file.dir , include.dirs = FALSE)
nexrad.L3.HCA.1.8deg.listfiles   <- list.files(path = nexrad.L3.HCA.1.8deg.file.dir , include.dirs = FALSE)
nexrad.L3.HCA.2.4deg.listfiles   <- list.files(path = nexrad.L3.HCA.2.4deg.file.dir , include.dirs = FALSE)
nexrad.L3.HCA.3.4deg.listfiles   <- list.files(path = nexrad.L3.HCA.3.4deg.file.dir , include.dirs = FALSE)

# define NEXRAD L3 SSTRUC data file path and names
#   file path
nexrad.L3.SSTRUC.file.dir        <- file.path(paste("/d1/serke/projects/HAIL_ACCUM_detect_ATEC/data/NEXRAD_L3/", case.yyyymmdd, "_", NEXRAD.radar, "/csv/alltimes_percase_files/", sep=""))
#   file name
#nexrad.L3.SSTRUC.file.name       <- "20120411_alltimes_KAMA_L3SSTRUC.csv"
nexrad.L3.SSTRUC.listfiles       <- list.files(path = nexrad.L3.SSTRUC.file.dir, include.dirs = FALSE)

#----------------------------------
# print some information to the GUI
#----------------------------------
if (flag.display.to.gui == 1) {
  print(paste("------------------------------------------------------------",                                                                        sep=""))
  print(paste("Case date  = ", case.yyyymmdd,                                                                                                        sep=""))
  print(paste("Case radar = ", NEXRAD.radar,                                                                                                         sep=""))
  print(paste("  Radar L2 data:",                                                                                                                    sep=""))
  print(paste("    There are ", length(nexrad.L2.listfiles), " nc format files in ", count.volumes, " radar volumes available for this case",        sep=""))
  print(paste("  Radar L3 HCA data:",                                                                                                                sep=""))
  print(paste("    There are ", length(nexrad.L3.HCA.0.5deg.listfiles), " nc format 0.5 deg files available for this case",                          sep=""))
  print(paste("    There are ", length(nexrad.L3.HCA.0.9deg.listfiles), " nc format 0.9 deg files available for this case",                          sep=""))
  print(paste("    There are ", length(nexrad.L3.HCA.1.5deg.listfiles), " nc format 1.5 deg files available for this case",                          sep=""))
  print(paste("    There are ", length(nexrad.L3.HCA.1.8deg.listfiles), " nc format 1.8 deg files available for this case",                          sep=""))
  print(paste("    There are ", length(nexrad.L3.HCA.2.4deg.listfiles), " nc format 2.4 deg files available for this case",                          sep=""))
  print(paste("    There are ", length(nexrad.L3.HCA.3.4deg.listfiles), " nc format 3.4 deg files available for this case",                          sep=""))
  print(paste("  Radar L3 Storm Tracking data:",                                                                                                     sep=""))
  print(paste("    There are ", length(nexrad.L3.SSTRUC.listfiles), " csv format files available for this case",                                     sep=""))
  print(paste("  Sounding data:",                                                                                                                    sep=""))
  print(paste("    There are ", length(balloon.sounding.listfiles), " csv format files available for this case",                                     sep=""))
  print(paste("    The times of these files are ", substr(balloon.sounding.file.name, 10, 11), ":", substr(balloon.sounding.file.name, 12, 13), "Z", sep=""))
} else {
  print("Display to screen flag is OFF")
}
#-------------------------------
# define NEXRAD L2 and L3 range arrays
#-------------------------------
#   range gate resolution 250 m
#   far field begins at 1000 m
nexrad.L2.range.array.m                <- seq(from = 1000, to = 458750, by = 250)   # [m]
nexrad.L3.range.array.m                <- seq(from = 1000, to = 300750, by = 250)   # [m], this only works if using the non-RADX version of L3 nc files

# define NEXRAD L2 azimuth array
nexrad.L2.azimuth.array.deg            <- seq(from = 0.00, to = 359.50, by = 0.50)  # [deg]

#-----------------------------
# ingest all input data sources
#-----------------------------
#   load the NEXRAD site location text file
NEXRAD.site.df                         <- read.csv(paste(nexrad.site.dataframetxt.dir, "nexrad_site.csv", sep = ""), header = FALSE, sep = ",", dec = ".", stringsAsFactors=FALSE)
colnames(NEXRAD.site.df)               <- c("NCDCID", "ICAO", "WBAN", "radname", "COUNTRY", "STATE", "COUNTY", "lat", "lon", "elev", "GMTdiff", "STN_TYPE")
head(NEXRAD.site.df)
#   load info specific to KAMA
NEXRAD.radar.wspace                    <- paste(" ", NEXRAD.radar, sep="")
ind.radar                              <- which(NEXRAD.site.df$ICAO == NEXRAD.radar.wspace)
radar.elev.ft                          <- NEXRAD.site.df$elev[ind.radar]
radar.elev.m                           <- NISTftTOmeter(radar.elev.ft)
radar.lat.deg                          <- NEXRAD.site.df$lat[ind.radar]
radar.lon.deg                          <- NEXRAD.site.df$lon[ind.radar]

#   ingest temp profile (closest/most recent sounding or model, 5 hours earlier at KAMA)
balloon.sounding.df                    <- read.csv(paste(balloon.sounding.file.dir, balloon.sounding.file.name, sep = ""), header = FALSE, sep = " ", dec = ".", stringsAsFactors=FALSE)
colnames(balloon.sounding.df)          <- c("Pres.hPa", "hgt.m", "temp.amb.C", "temp.dew.C" )
head(balloon.sounding.df)

#   ingest L3 storm structure (SSTRUC) csv files
# USE deparse(nexrad.L3.SSTRUC.alltimes.df$time.yyyymmddhhmmss.utc[1])) to see full yyyymmddhhmmss
nexrad.L3.SSTRUC.alltimes.df           <- read.csv(paste(nexrad.L3.SSTRUC.file.dir, nexrad.L3.SSTRUC.listfiles, sep=""), header = TRUE, sep = ",", dec = ".", stringsAsFactors=FALSE)
colnames(nexrad.L3.SSTRUC.alltimes.df) <- c("time.yyyymmddhhmmss.utc", "lon.deg", "lat.deg", "wsr.id", "cell.id", "range.nmi", "azim.deg", "hgt.base.kft", "hgt.top.kft", "vil.kgm2", "z.max.dbz", "hgt.z.max.kft" )
head(nexrad.L3.SSTRUC.alltimes.df)

# initialize polar coordinate sfc.hail.tot.accum.cm 1200 x 360 matrix
#   will be added to with each volume 'a'
sfc.hail.totaccum.cm                   <- array(data = 0L, dim = c(1200, 360))

# loop over 'a' number of volumes and load 'b' number of tilts for the given volume 
#for (a in 1:count.volumes) { 
for (a in 15:15) {
#a <- 11

  sfc.hail.1volaccum.cm          <- array(data = 0L, dim = c(1200, 360))
  
  for (b in count.files.in.volume.array[a]) {

    print(paste("--------------------------------------------------------------", sep=""))
    print(paste("Looping through L2 data volume #", a, " out of ", count.volumes, sep=""))
    
    if (flag.display.to.gui == 1) {
      print(paste("  VCP = ", nexrad.L2.tilt1.VCP, sep=""))
      print(paste("  Loading ", paste(nexrad.L2.file.dir, nexrad.L2.listfiles[1+a*b+(b-(b-0))], sep="")))
      print(paste("  Loading ", paste(nexrad.L2.file.dir, nexrad.L2.listfiles[1+a*b+(b-(b-1))], sep="")))
      print(paste("  Loading ", paste(nexrad.L2.file.dir, nexrad.L2.listfiles[1+a*b+(b-(b-2))], sep="")))
      print(paste("  Loading ", paste(nexrad.L2.file.dir, nexrad.L2.listfiles[1+a*b+(b-(b-3))], sep="")))
      print(paste("  Loading ", paste(nexrad.L2.file.dir, nexrad.L2.listfiles[1+a*b+(b-(b-4))], sep="")))
      print(paste("  Loading ", paste(nexrad.L2.file.dir, nexrad.L2.listfiles[1+a*b+(b-(b-5))], sep="")))
      print(paste("  Loading ", paste(nexrad.L2.file.dir, nexrad.L2.listfiles[1+a*b+(b-(b-6))], sep="")))
      print(paste("  Loading ", paste(nexrad.L2.file.dir, nexrad.L2.listfiles[1+a*b+(b-(b-7))], sep="")))
      print(paste("  Loading ", paste(nexrad.L2.file.dir, nexrad.L2.listfiles[1+a*b+(b-(b-8))], sep="")))
      print(paste("  Loading ", paste(nexrad.L2.file.dir, nexrad.L2.listfiles[1+a*b+(b-(b-9))], sep="")))
      print(paste("  Loading ", paste(nexrad.L2.file.dir, nexrad.L2.listfiles[1+a*b+(b-(b-10))], sep="")))
      print(paste("  Loading ", paste(nexrad.L2.file.dir, nexrad.L2.listfiles[1+a*b+(b-(b-11))], sep="")))
      print(paste("  Loading ", paste(nexrad.L2.file.dir, nexrad.L2.listfiles[1+a*b+(b-(b-12))], sep="")))
      print(paste("  Loading ", paste(nexrad.L2.file.dir, nexrad.L2.listfiles[1+a*b+(b-(b-13))], sep="")))
    } else {
      print("")
    }
    
    # calc tot.min = hh*60 + mm, use tilt3 as the representative time for this radar volume (since there are only the lowest 6 tilts for NEXRAD L3 files, dealt with later)
    tot.min.L2.tilt3 <- as.numeric(substr(nexrad.L2.listfiles[a*b+(b-(b-2))], 21, 22)) * 60 + as.numeric(substr(nexrad.L2.listfiles[1+a*b+(b-(b-2))], 23, 24))
    
    #   ingest NEXRAD L2 dBZ field data
    # NC TILT FILES ARE LOCATED AT: gturb3:/d1/serke/data/radarRaw/KAMA/20120411, about 12 tilts per volume x 5 volumes from 05:00 to 06:00 UTC
    nexrad.L2.tilt1.data          <- nc_open(paste(nexrad.L2.file.dir, nexrad.L2.listfiles[1+a*b+(b-(b-0))], sep=""), write = FALSE, verbose = FALSE) 
    #print(paste("The file has", nexrad.L2.tilt1.data$nvars, "variables"))
    nexrad.L2.tilt1.data.var.num  <- seq(1, nexrad.L2.tilt1.data$nvars, by=1)
    for (i in 1:length(nexrad.L2.tilt1.data.var.num)) {
      nexrad.L2.tilt1.data.nam    <- paste("v", nexrad.L2.tilt1.data.var.num[i], sep = "")
      assign(nexrad.L2.tilt1.data.nam, nexrad.L2.tilt1.data$var[[nexrad.L2.tilt1.data.var.num[i]]])
    }  # end of for (i in 1:length())
    if (nexrad.L2.tilt1.data$nvars == 44) {
      nexrad.L2.tilt1.AZIM  <- ncvar_get( nexrad.L2.tilt1.data, v38 )
      nexrad.L2.tilt1.Z     <- ncvar_get( nexrad.L2.tilt1.data, v42 )
    } else if (nexrad.L2.tilt1.data$nvars < 44) {
      print("Error")
    }  # end of if (data.visib ...)
    # also get VCP from data frame
    nexrad.L2.tilt1.VCP                 <- ncatt_get(nexrad.L2.tilt1.data, 0, "VCP")$value
    nc_close(nexrad.L2.tilt1.data)

    nexrad.L2.tilt2.data          <- nc_open(paste(nexrad.L2.file.dir, nexrad.L2.listfiles[1+a*b+(b-(b-1))], sep=""), write = FALSE, verbose = FALSE)
    #print(paste("The file has", nexrad.L2.tilt2.data$nvars, "variables"))
    nexrad.L2.tilt2.data.var.num  <- seq(1, nexrad.L2.tilt2.data$nvars, by=1)
    for (i in 1:length(nexrad.L2.tilt2.data.var.num)) {
      nexrad.L2.tilt2.data.nam    <- paste("v", nexrad.L2.tilt2.data.var.num[i], sep = "")
      assign(nexrad.L2.tilt2.data.nam, nexrad.L2.tilt2.data$var[[nexrad.L2.tilt2.data.var.num[i]]])
    }  # end of for (i in 1:length())
    if (nexrad.L2.tilt2.data$nvars == 44) {
      nexrad.L2.tilt2.AZIM  <- ncvar_get( nexrad.L2.tilt2.data, v38 )
      nexrad.L2.tilt2.Z     <- ncvar_get( nexrad.L2.tilt2.data, v42 )
    } else if (nexrad.L2.tilt2.data$nvars < 44) {
      print("Error")
    }  # end of if (data.visib ...)
    nc_close(nexrad.L2.tilt2.data)

    nexrad.L2.tilt3.data          <- nc_open(paste(nexrad.L2.file.dir, nexrad.L2.listfiles[1+a*b+(b-(b-2))], sep=""), write = FALSE, verbose = FALSE)
    #print(paste("The file has", nexrad.L2.tilt3.data$nvars, "variables"))
    nexrad.L2.tilt3.data.var.num  <- seq(1, nexrad.L2.tilt3.data$nvars, by=1)
    for (i in 1:length(nexrad.L2.tilt3.data.var.num)) {
      nexrad.L2.tilt3.data.nam    <- paste("v", nexrad.L2.tilt3.data.var.num[i], sep = "")
      assign(nexrad.L2.tilt3.data.nam, nexrad.L2.tilt3.data$var[[nexrad.L2.tilt3.data.var.num[i]]])
    }  # end of for (i in 1:length())
    if (nexrad.L2.tilt3.data$nvars == 44) {
      nexrad.L2.tilt3.AZIM  <- ncvar_get( nexrad.L2.tilt3.data, v38 )
      nexrad.L2.tilt3.Z     <- ncvar_get( nexrad.L2.tilt3.data, v42 )
    } else if (nexrad.L2.tilt3.data$nvars < 44) {
      print("Error")
    }  # end of if (data.visib ...)
    nc_close(nexrad.L2.tilt3.data)

    nexrad.L2.tilt4.data          <- nc_open(paste(nexrad.L2.file.dir, nexrad.L2.listfiles[1+a*b+(b-(b-3))], sep=""), write = FALSE, verbose = FALSE)
    #print(paste("The file has", nexrad.L2.tilt4.data$nvars, "variables"))
    nexrad.L2.tilt4.data.var.num  <- seq(1, nexrad.L2.tilt4.data$nvars, by=1)
    for (i in 1:length(nexrad.L2.tilt4.data.var.num)) {
      nexrad.L2.tilt4.data.nam    <- paste("v", nexrad.L2.tilt4.data.var.num[i], sep = "")
      assign(nexrad.L2.tilt4.data.nam, nexrad.L2.tilt4.data$var[[nexrad.L2.tilt4.data.var.num[i]]])
    }  # end of for (i in 1:length())
    if (nexrad.L2.tilt4.data$nvars == 44) {
      nexrad.L2.tilt4.AZIM  <- ncvar_get( nexrad.L2.tilt4.data, v38 )
      nexrad.L2.tilt4.Z     <- ncvar_get( nexrad.L2.tilt4.data, v42 )
    } else if (nexrad.L2.tilt4.data$nvars < 44) {
      print("Error")
    }  # end of if (data.visib ...)
    nc_close(nexrad.L2.tilt4.data)

    nexrad.L2.tilt5.data          <- nc_open(paste(nexrad.L2.file.dir, nexrad.L2.listfiles[1+a*b+(b-(b-4))], sep=""), write = FALSE, verbose = FALSE)
    #print(paste("The file has", nexrad.L2.tilt5.data$nvars, "variables"))
    nexrad.L2.tilt5.data.var.num  <- seq(1, nexrad.L2.tilt5.data$nvars, by=1)
    for (i in 1:length(nexrad.L2.tilt5.data.var.num)) {
      nexrad.L2.tilt5.data.nam    <- paste("v", nexrad.L2.tilt5.data.var.num[i], sep = "")
      assign(nexrad.L2.tilt5.data.nam, nexrad.L2.tilt5.data$var[[nexrad.L2.tilt5.data.var.num[i]]])
    }  # end of for (i in 1:length())
    if (nexrad.L2.tilt5.data$nvars == 44) {
      nexrad.L2.tilt5.AZIM  <- ncvar_get( nexrad.L2.tilt5.data, v38 )
      nexrad.L2.tilt5.Z     <- ncvar_get( nexrad.L2.tilt5.data, v42 )
    } else if (nexrad.L2.tilt5.data$nvars < 44) {
      print("Error")
    }  # end of if (data.visib ...)
    nc_close(nexrad.L2.tilt5.data)

    nexrad.L2.tilt6.data          <- nc_open(paste(nexrad.L2.file.dir, nexrad.L2.listfiles[1+a*b+(b-(b-5))], sep=""), write = FALSE, verbose = FALSE)
    #print(paste("The file has", nexrad.L2.tilt6.data$nvars, "variables"))
    nexrad.L2.tilt6.data.var.num  <- seq(1, nexrad.L2.tilt6.data$nvars, by=1)
    for (i in 1:length(nexrad.L2.tilt6.data.var.num)) {
      nexrad.L2.tilt6.data.nam    <- paste("v", nexrad.L2.tilt6.data.var.num[i], sep = "")
      assign(nexrad.L2.tilt6.data.nam, nexrad.L2.tilt6.data$var[[nexrad.L2.tilt6.data.var.num[i]]])
    }  # end of for (i in 1:length())
    if (nexrad.L2.tilt6.data$nvars == 44) {
      nexrad.L2.tilt6.AZIM  <- ncvar_get( nexrad.L2.tilt6.data, v38 )
      nexrad.L2.tilt6.Z     <- ncvar_get( nexrad.L2.tilt6.data, v42 )
    } else if (nexrad.L2.tilt6.data$nvars < 44) {
      print("Error")
    }  # end of if (data.visib ...)
    nc_close(nexrad.L2.tilt6.data)

    nexrad.L2.tilt7.data          <- nc_open(paste(nexrad.L2.file.dir, nexrad.L2.listfiles[1+a*b+(b-(b-6))], sep=""), write = FALSE, verbose = FALSE)
    #print(paste("The file has", nexrad.L2.tilt7.data$nvars, "variables"))
    nexrad.L2.tilt7.data.var.num  <- seq(1, nexrad.L2.tilt7.data$nvars, by=1)
    for (i in 1:length(nexrad.L2.tilt7.data.var.num)) {
      nexrad.L2.tilt7.data.nam    <- paste("v", nexrad.L2.tilt7.data.var.num[i], sep = "")
      assign(nexrad.L2.tilt7.data.nam, nexrad.L2.tilt7.data$var[[nexrad.L2.tilt7.data.var.num[i]]])
    }  # end of for (i in 1:length())
    if (nexrad.L2.tilt7.data$nvars == 44) {
      nexrad.L2.tilt7.AZIM  <- ncvar_get( nexrad.L2.tilt7.data, v38 )
      nexrad.L2.tilt7.Z     <- ncvar_get( nexrad.L2.tilt7.data, v42 )
    } else if (nexrad.L2.tilt7.data$nvars < 44) {
      print("Error")
    }  # end of if (data.visib ...)
    nc_close(nexrad.L2.tilt7.data)

    nexrad.L2.tilt8.data          <- nc_open(paste(nexrad.L2.file.dir, nexrad.L2.listfiles[1+a*b+(b-(b-7))], sep=""), write = FALSE, verbose = FALSE)
    #print(paste("The file has", nexrad.L2.tilt8.data$nvars, "variables"))
    nexrad.L2.tilt8.data.var.num  <- seq(1, nexrad.L2.tilt8.data$nvars, by=1)
    for (i in 1:length(nexrad.L2.tilt8.data.var.num)) {
      nexrad.L2.tilt8.data.nam    <- paste("v", nexrad.L2.tilt8.data.var.num[i], sep = "")
      assign(nexrad.L2.tilt8.data.nam, nexrad.L2.tilt8.data$var[[nexrad.L2.tilt8.data.var.num[i]]])
    }  # end of for (i in 1:length())
    if (nexrad.L2.tilt8.data$nvars == 44) {
      nexrad.L2.tilt8.AZIM  <- ncvar_get( nexrad.L2.tilt8.data, v38 )
      nexrad.L2.tilt8.Z     <- ncvar_get( nexrad.L2.tilt8.data, v42 )
    } else if (nexrad.L2.tilt8.data$nvars < 44) {
      print("Error")
    }  # end of if (data.visib ...)
    nc_close(nexrad.L2.tilt8.data)

    nexrad.L2.tilt9.data          <- nc_open(paste(nexrad.L2.file.dir, nexrad.L2.listfiles[1+a*b+(b-(b-8))], sep=""), write = FALSE, verbose = FALSE)
    #print(paste("The file has", nexrad.L2.tilt9.data$nvars, "variables"))
    nexrad.L2.tilt9.data.var.num  <- seq(1, nexrad.L2.tilt9.data$nvars, by=1)
    for (i in 1:length(nexrad.L2.tilt9.data.var.num)) {
      nexrad.L2.tilt9.data.nam    <- paste("v", nexrad.L2.tilt9.data.var.num[i], sep = "")
      assign(nexrad.L2.tilt9.data.nam, nexrad.L2.tilt9.data$var[[nexrad.L2.tilt9.data.var.num[i]]])
    }  # end of for (i in 1:length())
    if (nexrad.L2.tilt9.data$nvars == 44) {
      nexrad.L2.tilt9.AZIM  <- ncvar_get( nexrad.L2.tilt9.data, v38 )
      nexrad.L2.tilt9.Z     <- ncvar_get( nexrad.L2.tilt9.data, v42 )
    } else if (nexrad.L2.tilt9.data$nvars < 44) {
      print("Error")
    }  # end of if (data.visib ...)
    nc_close(nexrad.L2.tilt9.data)

    nexrad.L2.tilt10.data          <- nc_open(paste(nexrad.L2.file.dir, nexrad.L2.listfiles[1+a*b+(b-(b-9))], sep=""), write = FALSE, verbose = FALSE)
    #print(paste("The file has", nexrad.L2.tilt10.data$nvars, "variables"))
    nexrad.L2.tilt10.data.var.num  <- seq(1, nexrad.L2.tilt10.data$nvars, by=1)
    for (i in 1:length(nexrad.L2.tilt10.data.var.num)) {
      nexrad.L2.tilt10.data.nam    <- paste("v", nexrad.L2.tilt10.data.var.num[i], sep = "")
      assign(nexrad.L2.tilt10.data.nam, nexrad.L2.tilt10.data$var[[nexrad.L2.tilt10.data.var.num[i]]])
    }  # end of for (i in 1:length())
    if (nexrad.L2.tilt10.data$nvars == 44) {
      nexrad.L2.tilt10.AZIM  <- ncvar_get( nexrad.L2.tilt10.data, v38 )
      nexrad.L2.tilt10.Z     <- ncvar_get( nexrad.L2.tilt10.data, v42 )
    } else if (nexrad.L2.tilt10.data$nvars < 44) {
      print("Error")
    }  # end of if (data.visib ...)
    nc_close(nexrad.L2.tilt10.data)

    nexrad.L2.tilt11.data          <- nc_open(paste(nexrad.L2.file.dir, nexrad.L2.listfiles[1+a*b+(b-(b-10))], sep=""), write = FALSE, verbose = FALSE)
    #print(paste("The file has", nexrad.L2.tilt11.data$nvars, "variables"))
    nexrad.L2.tilt11.data.var.num  <- seq(1, nexrad.L2.tilt11.data$nvars, by=1)
    for (i in 1:length(nexrad.L2.tilt11.data.var.num)) {
      nexrad.L2.tilt11.data.nam    <- paste("v", nexrad.L2.tilt11.data.var.num[i], sep = "")
      assign(nexrad.L2.tilt11.data.nam, nexrad.L2.tilt11.data$var[[nexrad.L2.tilt11.data.var.num[i]]])
    }  # end of for (i in 1:length())
    if (nexrad.L2.tilt11.data$nvars == 44) {
      nexrad.L2.tilt11.AZIM  <- ncvar_get( nexrad.L2.tilt11.data, v38 )
      nexrad.L2.tilt11.Z     <- ncvar_get( nexrad.L2.tilt11.data, v42 )
    } else if (nexrad.L2.tilt11.data$nvars < 44) {
      print("Error")
    }  # end of if (data.visib ...)
    nc_close(nexrad.L2.tilt11.data)

    nexrad.L2.tilt12.data          <- nc_open(paste(nexrad.L2.file.dir, nexrad.L2.listfiles[1+a*b+(b-(b-11))], sep=""), write = FALSE, verbose = FALSE)
    #print(paste("The file has", nexrad.L2.tilt12.data$nvars, "variables"))
    nexrad.L2.tilt12.data.var.num  <- seq(1, nexrad.L2.tilt12.data$nvars, by=1)
    for (i in 1:length(nexrad.L2.tilt12.data.var.num)) {
      nexrad.L2.tilt12.data.nam <- paste("v", nexrad.L2.tilt12.data.var.num[i], sep = "")
      assign(nexrad.L2.tilt12.data.nam, nexrad.L2.tilt12.data$var[[nexrad.L2.tilt12.data.var.num[i]]])
    }  # end of for (i in 1:length())
    if (nexrad.L2.tilt12.data$nvars == 44) {
      nexrad.L2.tilt12.AZIM     <- ncvar_get( nexrad.L2.tilt12.data, v38 )
      nexrad.L2.tilt12.Z        <- ncvar_get( nexrad.L2.tilt12.data, v42 )
    } else if (nexrad.L2.tilt12.data$nvars < 44) {
      print("Error")
    }  # end of if (data.visib ...)
    nc_close(nexrad.L2.tilt12.data)

    nexrad.L2.tilt13.data          <- nc_open(paste(nexrad.L2.file.dir, nexrad.L2.listfiles[1+a*b+(b-(b-12))], sep=""), write = FALSE, verbose = FALSE)
    #print(paste("The file has", nexrad.L2.tilt13.data$nvars, "variables"))
    nexrad.L2.tilt13.data.var.num  <- seq(1, nexrad.L2.tilt13.data$nvars, by=1)
    for (i in 1:length(nexrad.L2.tilt13.data.var.num)) {
      nexrad.L2.tilt13.data.nam <- paste("v", nexrad.L2.tilt13.data.var.num[i], sep = "")
      assign(nexrad.L2.tilt13.data.nam, nexrad.L2.tilt13.data$var[[nexrad.L2.tilt13.data.var.num[i]]])
    }  # end of for (i in 1:length())
    if (nexrad.L2.tilt13.data$nvars == 44) {
      nexrad.L2.tilt13.AZIM     <- ncvar_get( nexrad.L2.tilt13.data, v38 )
      nexrad.L2.tilt13.Z        <- ncvar_get( nexrad.L2.tilt13.data, v42 )
    } else if (nexrad.L2.tilt13.data$nvars < 44) {
      print("Error")
    }  # end of if (data.visib ...)
    nc_close(nexrad.L2.tilt13.data)

    nexrad.L2.tilt14.data          <- nc_open(paste(nexrad.L2.file.dir, nexrad.L2.listfiles[1+a*b+(b-(b-13))], sep=""), write = FALSE, verbose = FALSE)
    #print(paste("The file has", nexrad.L2.tilt14.data$nvars, "variables"))
    nexrad.L2.tilt14.data.var.num  <- seq(1, nexrad.L2.tilt14.data$nvars, by=1)
    for (i in 1:length(nexrad.L2.tilt14.data.var.num)) {
      nexrad.L2.tilt14.data.nam <- paste("v", nexrad.L2.tilt14.data.var.num[i], sep = "")
      assign(nexrad.L2.tilt14.data.nam, nexrad.L2.tilt14.data$var[[nexrad.L2.tilt14.data.var.num[i]]])
    }  # end of for (i in 1:length())
    if (nexrad.L2.tilt14.data$nvars == 44) {
      nexrad.L2.tilt14.AZIM     <- ncvar_get( nexrad.L2.tilt14.data, v38 )
      nexrad.L2.tilt14.Z        <- ncvar_get( nexrad.L2.tilt14.data, v42 )
    } else if (nexrad.L2.tilt14.data$nvars < 44) {
      print("Error")
    }  # end of if (data.visib ...)
    nc_close(nexrad.L2.tilt14.data)
    
  }  # end of for (b in ...)

  # identify the closest L3 SSTRUC time, and then read in all rows with similar times
  print(paste("Identifying L3 SSTRUC closest temporal file match", sep=""))
  tot.min.L3.SSTRUC.files     <- (as.numeric(substr(as.character(nexrad.L3.SSTRUC.alltimes.df$time.yyyymmddhhmmss.utc), 9, 10)) * 60) + as.numeric(substr(as.character(nexrad.L3.SSTRUC.alltimes.df$time.yyyymmddhhmmss.utc), 11, 12))
  ind.L3.SSTRUC.closest       <- which.min(abs(tot.min.L2.tilt3 - tot.min.L3.SSTRUC.files))
  ind.L3.SSTRUC.all.closest   <- which(nexrad.L3.SSTRUC.alltimes.df$time.yyyymmddhhmmss.utc == nexrad.L3.SSTRUC.alltimes.df$time.yyyymmddhhmmss.utc[ind.L3.SSTRUC.closest])
  
  # identify the closes L3 HCA file (in time) for each of the 6 tilt angles
  print(paste("Identifying L3 HCA closest temporal file match", sep=""))
  tot.min.L3.HCA.0.5deg.files <- (as.numeric(substr(nexrad.L3.HCA.0.5deg.listfiles, 10, 11)) * 60) + as.numeric(substr(nexrad.L3.HCA.0.5deg.listfiles, 12, 13))
  tot.min.L3.HCA.0.9deg.files <- (as.numeric(substr(nexrad.L3.HCA.0.9deg.listfiles, 10, 11)) * 60) + as.numeric(substr(nexrad.L3.HCA.0.9deg.listfiles, 12, 13))
  tot.min.L3.HCA.1.5deg.files <- (as.numeric(substr(nexrad.L3.HCA.1.5deg.listfiles, 10, 11)) * 60) + as.numeric(substr(nexrad.L3.HCA.1.5deg.listfiles, 12, 13))
  tot.min.L3.HCA.1.8deg.files <- (as.numeric(substr(nexrad.L3.HCA.1.8deg.listfiles, 10, 11)) * 60) + as.numeric(substr(nexrad.L3.HCA.1.8deg.listfiles, 12, 13))
  tot.min.L3.HCA.2.4deg.files <- (as.numeric(substr(nexrad.L3.HCA.2.4deg.listfiles, 10, 11)) * 60) + as.numeric(substr(nexrad.L3.HCA.2.4deg.listfiles, 12, 13))
  tot.min.L3.HCA.3.4deg.files <- (as.numeric(substr(nexrad.L3.HCA.3.4deg.listfiles, 10, 11)) * 60) + as.numeric(substr(nexrad.L3.HCA.3.4deg.listfiles, 12, 13))
  ind.L3.HCA.0.5deg.closest   <- which.min(abs(tot.min.L2.tilt3 - tot.min.L3.HCA.0.5deg.files))
  ind.L3.HCA.0.9deg.closest   <- which.min(abs(tot.min.L2.tilt3 - tot.min.L3.HCA.0.9deg.files))
  ind.L3.HCA.1.5deg.closest   <- which.min(abs(tot.min.L2.tilt3 - tot.min.L3.HCA.1.5deg.files))
  ind.L3.HCA.1.8deg.closest   <- which.min(abs(tot.min.L2.tilt3 - tot.min.L3.HCA.1.8deg.files))
  ind.L3.HCA.2.4deg.closest   <- which.min(abs(tot.min.L2.tilt3 - tot.min.L3.HCA.2.4deg.files))
  ind.L3.HCA.3.4deg.closest   <- which.min(abs(tot.min.L2.tilt3 - tot.min.L3.HCA.3.4deg.files))
  
  # ingest NEXRAD L3 HCA (Park et al, 2009) and use Hail with rain category (cat = 9 and 10)
  #   for tilt 1
  nexrad.L3.tilt1.data          <- nc_open(paste(nexrad.L3.HCA.0.5deg.file.dir, nexrad.L3.HCA.0.5deg.listfiles[ind.L3.HCA.0.5deg.closest], sep=""), write = FALSE, verbose = FALSE)
  print(paste("  Loading ", paste(nexrad.L3.HCA.0.5deg.file.dir, nexrad.L3.HCA.0.5deg.listfiles[ind.L3.HCA.0.5deg.closest], sep="")))
  #print(paste("The file has", nexrad.L3.tilt1.data$nvars, "variables"))
  nexrad.L3.tilt1.data.var.num  <- seq(1, nexrad.L3.tilt1.data$nvars, by=1)
  for (i in 1:length(nexrad.L3.tilt1.data.var.num)) {
    nexrad.L3.tilt1.data.nam    <- paste("v", nexrad.L3.tilt1.data.var.num[i], sep = "")
    assign(nexrad.L3.tilt1.data.nam, nexrad.L3.tilt1.data$var[[nexrad.L3.tilt1.data.var.num[i]]])
  }  # end of for (i in 1:length())
  if (nexrad.L3.tilt1.data$nvars == 6) {
    nexrad.L3.tilt1.LAT  <- ncvar_get( nexrad.L3.tilt1.data, v3 )
    nexrad.L3.tilt1.LON  <- ncvar_get( nexrad.L3.tilt1.data, v4 )
    nexrad.L3.tilt1.HCA  <- ncvar_get( nexrad.L3.tilt1.data, v6 )
  } else if (nexrad.L3.tilt1.data$nvars < 6) {
    print("Error")
  }  # end of if (data.visib ...)
  nc_close(nexrad.L3.tilt1.data)

  #   for tilt 2
  nexrad.L3.tilt2.data          <- nc_open(paste(nexrad.L3.HCA.0.9deg.file.dir, nexrad.L3.HCA.0.9deg.listfiles[ind.L3.HCA.0.9deg.closest], sep=""), write = FALSE, verbose = FALSE)
  print(paste("  Loading ", paste(nexrad.L3.HCA.0.9deg.file.dir, nexrad.L3.HCA.0.9deg.listfiles[ind.L3.HCA.0.9deg.closest], sep="")))
  #print(paste("The file has", nexrad.L3.tilt2.data$nvars, "variables"))
  nexrad.L3.tilt2.data.var.num  <- seq(1, nexrad.L3.tilt2.data$nvars, by=1)
  for (i in 1:length(nexrad.L3.tilt2.data.var.num)) {
    nexrad.L3.tilt2.data.nam    <- paste("v", nexrad.L3.tilt2.data.var.num[i], sep = "")
    assign(nexrad.L3.tilt2.data.nam, nexrad.L3.tilt2.data$var[[nexrad.L3.tilt2.data.var.num[i]]])
  }  # end of for (i in 1:length())
  if (nexrad.L3.tilt2.data$nvars == 6) {
    nexrad.L3.tilt2.LAT  <- ncvar_get( nexrad.L3.tilt2.data, v3 )
    nexrad.L3.tilt2.LON  <- ncvar_get( nexrad.L3.tilt2.data, v4 )
    nexrad.L3.tilt2.HCA  <- ncvar_get( nexrad.L3.tilt2.data, v6 )
  } else if (nexrad.L3.tilt2.data$nvars < 6) {
    print("Error")
  }  # end of if (data.visib ...)
  nc_close(nexrad.L3.tilt2.data)

  #   for tilt 3
  nexrad.L3.tilt3.data          <- nc_open(paste(nexrad.L3.HCA.1.5deg.file.dir, nexrad.L3.HCA.1.5deg.listfiles[ind.L3.HCA.1.5deg.closest], sep=""), write = FALSE, verbose = FALSE)
  print(paste("  Loading ", paste(nexrad.L3.HCA.1.5deg.file.dir, nexrad.L3.HCA.1.5deg.listfiles[ind.L3.HCA.1.5deg.closest], sep="")))
  #print(paste("The file has", nexrad.L3.tilt3.data$nvars, "variables"))
  nexrad.L3.tilt3.data.var.num  <- seq(1, nexrad.L3.tilt3.data$nvars, by=1)
  for (i in 1:length(nexrad.L3.tilt3.data.var.num)) {
    nexrad.L3.tilt3.data.nam    <- paste("v", nexrad.L3.tilt3.data.var.num[i], sep = "")
    assign(nexrad.L3.tilt3.data.nam, nexrad.L3.tilt3.data$var[[nexrad.L3.tilt3.data.var.num[i]]])
  }  # end of for (i in 1:length())
  if (nexrad.L3.tilt3.data$nvars == 6) {
    nexrad.L3.tilt3.LAT  <- ncvar_get( nexrad.L3.tilt3.data, v3 )
    nexrad.L3.tilt3.LON  <- ncvar_get( nexrad.L3.tilt3.data, v4 )
    nexrad.L3.tilt3.HCA  <- ncvar_get( nexrad.L3.tilt3.data, v6 )
  } else if (nexrad.L3.tilt3.data$nvars < 6) {
    print("Error")
  }  # end of if (data.visib ...)
  nc_close(nexrad.L3.tilt3.data)

  #   for tilt 4
  nexrad.L3.tilt4.data          <- nc_open(paste(nexrad.L3.HCA.1.8deg.file.dir, nexrad.L3.HCA.1.8deg.listfiles[ind.L3.HCA.1.8deg.closest], sep=""), write = FALSE, verbose = FALSE)
  print(paste("  Loading ", paste(nexrad.L3.HCA.1.8deg.file.dir, nexrad.L3.HCA.1.8deg.listfiles[ind.L3.HCA.1.8deg.closest], sep="")))
  #print(paste("The file has", nexrad.L3.tilt4.data$nvars, "variables"))
  nexrad.L3.tilt4.data.var.num  <- seq(1, nexrad.L3.tilt4.data$nvars, by=1)
  for (i in 1:length(nexrad.L3.tilt4.data.var.num)) {
    nexrad.L3.tilt4.data.nam    <- paste("v", nexrad.L3.tilt4.data.var.num[i], sep = "")
    assign(nexrad.L3.tilt4.data.nam, nexrad.L3.tilt4.data$var[[nexrad.L3.tilt4.data.var.num[i]]])
  }  # end of for (i in 1:length())
  if (nexrad.L3.tilt4.data$nvars == 6) {
    nexrad.L3.tilt4.LAT  <- ncvar_get( nexrad.L3.tilt4.data, v3 )
    nexrad.L3.tilt4.LON  <- ncvar_get( nexrad.L3.tilt4.data, v4 )
    nexrad.L3.tilt4.HCA  <- ncvar_get( nexrad.L3.tilt4.data, v6 )
  } else if (nexrad.L3.tilt4.data$nvars < 6) {
    print("Error")
  }  # end of if (data.visib ...)
  nc_close(nexrad.L3.tilt4.data)

  #   for tilt 5
  nexrad.L3.tilt5.data          <- nc_open(paste(nexrad.L3.HCA.2.4deg.file.dir, nexrad.L3.HCA.2.4deg.listfiles[ind.L3.HCA.2.4deg.closest], sep=""), write = FALSE, verbose = FALSE)
  print(paste("  Loading ", paste(nexrad.L3.HCA.2.4deg.file.dir, nexrad.L3.HCA.2.4deg.listfiles[ind.L3.HCA.2.4deg.closest], sep="")))
  #print(paste("The file has", nexrad.L3.tilt5.data$nvars, "variables"))
  nexrad.L3.tilt5.data.var.num  <- seq(1, nexrad.L3.tilt5.data$nvars, by=1)
  for (i in 1:length(nexrad.L3.tilt5.data.var.num)) {
    nexrad.L3.tilt5.data.nam    <- paste("v", nexrad.L3.tilt5.data.var.num[i], sep = "")
    assign(nexrad.L3.tilt5.data.nam, nexrad.L3.tilt5.data$var[[nexrad.L3.tilt5.data.var.num[i]]])
  }  # end of for (i in 1:length())
  if (nexrad.L3.tilt5.data$nvars == 6) {
    nexrad.L3.tilt5.LAT  <- ncvar_get( nexrad.L3.tilt5.data, v3 )
    nexrad.L3.tilt5.LON  <- ncvar_get( nexrad.L3.tilt5.data, v4 )
    nexrad.L3.tilt5.HCA  <- ncvar_get( nexrad.L3.tilt5.data, v6 )
  } else if (nexrad.L3.tilt5.data$nvars < 6) {
    print("Error")
  }  # end of if (data.visib ...)
  nc_close(nexrad.L3.tilt5.data)

  #   for tilt 6
  nexrad.L3.tilt6.data          <- nc_open(paste(nexrad.L3.HCA.3.4deg.file.dir, nexrad.L3.HCA.3.4deg.listfiles[ind.L3.HCA.3.4deg.closest], sep=""), write = FALSE, verbose = FALSE)
  print(paste("  Loading ", paste(nexrad.L3.HCA.2.4deg.file.dir, nexrad.L3.HCA.3.4deg.listfiles[ind.L3.HCA.2.4deg.closest], sep="")))
  #print(paste("The file has", nexrad.L3.tilt6.data$nvars, "variables"))
  nexrad.L3.tilt6.data.var.num  <- seq(1, nexrad.L3.tilt6.data$nvars, by=1)
  for (i in 1:length(nexrad.L3.tilt6.data.var.num)) {
    nexrad.L3.tilt6.data.nam    <- paste("v", nexrad.L3.tilt6.data.var.num[i], sep = "")
    assign(nexrad.L3.tilt6.data.nam, nexrad.L3.tilt6.data$var[[nexrad.L3.tilt6.data.var.num[i]]])
  }  # end of for (i in 1:length())
  if (nexrad.L3.tilt6.data$nvars == 6) {
    nexrad.L3.tilt6.LAT  <- ncvar_get( nexrad.L3.tilt6.data, v3 )
    nexrad.L3.tilt6.LON  <- ncvar_get( nexrad.L3.tilt6.data, v4 )
    nexrad.L3.tilt6.HCA  <- ncvar_get( nexrad.L3.tilt6.data, v6 )
  } else if (nexrad.L3.tilt6.data$nvars < 6) {
    print("Error")
  }  # end of if (data.visib ...)
  nc_close(nexrad.L3.tilt6.data)

  #----------------------------
  # Manipulate the ingested data
  #----------------------------
  print(paste("------------------------------------------", sep=""))
  print(paste("Manipulate ingested data streams as needed", sep=""))

  # NEXRAD VCP volume completion time (time.delta [m]) lookup table
  #   time.delta is defined as time between sucessive radar scans [min]
  print(paste("  Consult lookup table for VCP time.delta"))
  if (nexrad.L2.tilt1.VCP == 12) {
    time.delta.m <- 4.1
  } else if (nexrad.L2.tilt1.VCP == 212) {
    time.delta.m <- 4.5
  } else if (nexrad.L2.tilt1.VCP == 215 | nexrad.L2.tilt1.VCP == 121) {
    time.delta.m <- 6.0
  } else if (nexrad.L2.tilt1.VCP == 31 | nexrad.L2.tilt1.VCP == 32 | nexrad.L2.tilt1.VCP == 35) {
    time.delta.m <- 10.0
  } else {
    print("      Error: given VCP not accounted for")
  }
  print(paste("    Nexrad L2 volume time.delta.m based on VCP = ", nexrad.L2.tilt1.VCP, " is ", time.delta.m, " min"))
  
  # create new L2 Z matrices from raw input.  Apparently, Z scale folds above 30.5 to negative numbers
  print(paste("  Unfold L2 Z values above 30.5 dBZ", sep=""))
  nexrad.L2.tilt1.Z.new                   <- nexrad.L2.tilt1.Z
  nexrad.L2.tilt2.Z.new                   <- nexrad.L2.tilt2.Z
  nexrad.L2.tilt3.Z.new                   <- nexrad.L2.tilt3.Z
  nexrad.L2.tilt4.Z.new                   <- nexrad.L2.tilt4.Z
  nexrad.L2.tilt5.Z.new                   <- nexrad.L2.tilt5.Z
  nexrad.L2.tilt6.Z.new                   <- nexrad.L2.tilt6.Z
  nexrad.L2.tilt7.Z.new                   <- nexrad.L2.tilt7.Z
  nexrad.L2.tilt8.Z.new                   <- nexrad.L2.tilt8.Z
  nexrad.L2.tilt9.Z.new                   <- nexrad.L2.tilt9.Z
  nexrad.L2.tilt10.Z.new                  <- nexrad.L2.tilt10.Z
  nexrad.L2.tilt11.Z.new                  <- nexrad.L2.tilt11.Z
  nexrad.L2.tilt12.Z.new                  <- nexrad.L2.tilt12.Z
  nexrad.L2.tilt13.Z.new                  <- nexrad.L2.tilt13.Z
  nexrad.L2.tilt14.Z.new                  <- nexrad.L2.tilt14.Z
  ind.tilt1.lo.Z                          <- which(nexrad.L2.tilt1.Z  < -35)
  ind.tilt2.lo.Z                          <- which(nexrad.L2.tilt2.Z  < -35)
  ind.tilt3.lo.Z                          <- which(nexrad.L2.tilt3.Z  < -35)
  ind.tilt4.lo.Z                          <- which(nexrad.L2.tilt4.Z  < -35)
  ind.tilt5.lo.Z                          <- which(nexrad.L2.tilt5.Z  < -35)
  ind.tilt6.lo.Z                          <- which(nexrad.L2.tilt6.Z  < -35)
  ind.tilt7.lo.Z                          <- which(nexrad.L2.tilt7.Z  < -35)
  ind.tilt8.lo.Z                          <- which(nexrad.L2.tilt8.Z  < -35)
  ind.tilt9.lo.Z                          <- which(nexrad.L2.tilt9.Z  < -35)
  ind.tilt10.lo.Z                         <- which(nexrad.L2.tilt10.Z < -35)
  ind.tilt11.lo.Z                         <- which(nexrad.L2.tilt11.Z < -35)
  ind.tilt12.lo.Z                         <- which(nexrad.L2.tilt12.Z < -35)
  ind.tilt13.lo.Z                         <- which(nexrad.L2.tilt13.Z < -35)
  ind.tilt14.lo.Z                         <- which(nexrad.L2.tilt14.Z < -35)
  nexrad.L2.tilt1.Z.new[ind.tilt1.lo.Z]   <- nexrad.L2.tilt1.Z[ind.tilt1.lo.Z]   + 99 + 30.5
  nexrad.L2.tilt2.Z.new[ind.tilt2.lo.Z]   <- nexrad.L2.tilt2.Z[ind.tilt2.lo.Z]   + 99 + 30.5
  nexrad.L2.tilt3.Z.new[ind.tilt3.lo.Z]   <- nexrad.L2.tilt3.Z[ind.tilt3.lo.Z]   + 99 + 30.5
  nexrad.L2.tilt4.Z.new[ind.tilt4.lo.Z]   <- nexrad.L2.tilt4.Z[ind.tilt4.lo.Z]   + 99 + 30.5
  nexrad.L2.tilt5.Z.new[ind.tilt5.lo.Z]   <- nexrad.L2.tilt5.Z[ind.tilt5.lo.Z]   + 99 + 30.5
  nexrad.L2.tilt6.Z.new[ind.tilt6.lo.Z]   <- nexrad.L2.tilt6.Z[ind.tilt6.lo.Z]   + 99 + 30.5
  nexrad.L2.tilt7.Z.new[ind.tilt7.lo.Z]   <- nexrad.L2.tilt7.Z[ind.tilt7.lo.Z]   + 99 + 30.5
  nexrad.L2.tilt8.Z.new[ind.tilt8.lo.Z]   <- nexrad.L2.tilt8.Z[ind.tilt8.lo.Z]   + 99 + 30.5
  nexrad.L2.tilt9.Z.new[ind.tilt9.lo.Z]   <- nexrad.L2.tilt9.Z[ind.tilt9.lo.Z]   + 99 + 30.5
  nexrad.L2.tilt10.Z.new[ind.tilt10.lo.Z] <- nexrad.L2.tilt10.Z[ind.tilt10.lo.Z] + 99 + 30.5
  nexrad.L2.tilt11.Z.new[ind.tilt11.lo.Z] <- nexrad.L2.tilt11.Z[ind.tilt11.lo.Z] + 99 + 30.5
  nexrad.L2.tilt12.Z.new[ind.tilt12.lo.Z] <- nexrad.L2.tilt12.Z[ind.tilt12.lo.Z] + 99 + 30.5
  nexrad.L2.tilt13.Z.new[ind.tilt13.lo.Z] <- nexrad.L2.tilt13.Z[ind.tilt13.lo.Z] + 99 + 30.5
  nexrad.L2.tilt14.Z.new[ind.tilt14.lo.Z] <- nexrad.L2.tilt14.Z[ind.tilt14.lo.Z] + 99 + 30.5

  # load first 6 tilts of Z.new into a 3D matrix
  nexrad.L2.tilt1to6.Z.new                <- array(data = NA, dim = c(dim(nexrad.L2.tilt1.Z.new)[1], dim(nexrad.L2.tilt1.Z.new)[2], length(theta.radar.volume.L3.deg)))
  nexrad.L2.tilt1to6.Z.new[,,1]           <- nexrad.L2.tilt1.Z.new
  nexrad.L2.tilt2.Z.new                   <- rbind(nexrad.L2.tilt2.Z.new, array(data = NA, dim = c(dim(nexrad.L2.tilt1.Z.new)[1]-dim(nexrad.L2.tilt2.Z.new)[1], dim(nexrad.L2.tilt2.Z.new)[2])))
  nexrad.L2.tilt1to6.Z.new[,,2]           <- nexrad.L2.tilt2.Z.new
  nexrad.L2.tilt3.Z.new                   <- rbind(nexrad.L2.tilt3.Z.new, array(data = NA, dim = c(dim(nexrad.L2.tilt1.Z.new)[1]-dim(nexrad.L2.tilt3.Z.new)[1], dim(nexrad.L2.tilt3.Z.new)[2])))
  nexrad.L2.tilt1to6.Z.new[,,3]           <- nexrad.L2.tilt3.Z.new
  nexrad.L2.tilt4.Z.new                   <- rbind(nexrad.L2.tilt4.Z.new, array(data = NA, dim = c(dim(nexrad.L2.tilt1.Z.new)[1]-dim(nexrad.L2.tilt4.Z.new)[1], dim(nexrad.L2.tilt4.Z.new)[2])))
  nexrad.L2.tilt1to6.Z.new[,,4]           <- nexrad.L2.tilt4.Z.new
  nexrad.L2.tilt5.Z.new                   <- rbind(nexrad.L2.tilt5.Z.new, array(data = NA, dim = c(dim(nexrad.L2.tilt1.Z.new)[1]-dim(nexrad.L2.tilt5.Z.new)[1], dim(nexrad.L2.tilt5.Z.new)[2])))
  nexrad.L2.tilt1to6.Z.new[,,5]           <- nexrad.L2.tilt5.Z.new
  nexrad.L2.tilt6.Z.new                   <- rbind(nexrad.L2.tilt6.Z.new, array(data = NA, dim = c(dim(nexrad.L2.tilt1.Z.new)[1]-dim(nexrad.L2.tilt6.Z.new)[1], dim(nexrad.L2.tilt6.Z.new)[2])))
  nexrad.L2.tilt1to6.Z.new[,,6]           <- nexrad.L2.tilt6.Z.new

  #  calculate 3D matrix values for H.ARL.L2.m        
  #    is height above radar level from L2 data volume, in m
  #    accounts for height over curved surface
  #    matrix dim is [length(ranges), length(azimuths), length(tilts)]
  
  print(paste("  Calc L2 H.ARL, a 3D matrix", sep=""))
  H.ARL.L2.m                   <- array(data = 0L, dim = c(dim(nexrad.L2.tilt1.Z.new)[1], dim(nexrad.L2.tilt1.Z.new)[2], length(theta.radar.volume.L2.deg)))
  for (k in 1:length(theta.radar.volume.L2.deg)) {
    for (i in 1:length(nexrad.L2.range.array.m)) {
      
      H.ARL.L2.m[i,,k] <- nexrad.L2.range.array.m[i] * sin(NISTdegTOradian(theta.radar.volume.L2.deg[k])) + ((nexrad.L2.range.array.m[i]^2) / (2 * refract.index * (radius.earth.km * kmTOm)))
      
    }
  }    # end of for (k in 1:...)
  
  #  calculate matrix values for H.ARL.L3.m        
  #    is height above radar level from L3 data volume, in m
  #    accounts for height over curved surface
  #    matrix dim is [length(ranges), length(tilts)]
  print(paste("  Calc L3 H.ARL, a 3D matrix", sep=""))
  H.ARL.L3.m                      <- array(data = 0L, dim = c(dim(nexrad.L2.tilt1.Z.new)[1], dim(nexrad.L2.tilt1.Z.new)[2], length(theta.radar.volume.L3.deg)))
  for (k in 1:length(theta.radar.volume.L3.deg)) {
    for (i in 1:dim(nexrad.L2.tilt1.Z.new)[1]) {
      
      H.ARL.L3.m[i,,k] <- nexrad.L3.range.array.m[i] * sin(NISTdegTOradian(theta.radar.volume.L3.deg[k])) + ((nexrad.L3.range.array.m[i]^2) / (2 * refract.index * (radius.earth.km * kmTOm)))
      
    }
  }    # end of for (k in 1:...)
  
  # Calculate dist from radar along sfc (D.from.radar.m) 3D matrix
  print(paste("  Calc L3 horiz dist from radar, a 3D matrix"))
  Dist.from.radar.m                                <- array(data = NA, dim=c(dim(nexrad.L2.tilt1.Z.new)[1], dim(nexrad.L2.tilt1.Z.new)[2], length(theta.radar.volume.L3.deg)))
  for (k in 1:length(theta.radar.volume.L3.deg)) {
    for (i in 1:dim(nexrad.L2.tilt1.Z.new)[1]) {
      Dist.from.radar.m[i,,k] <- nexrad.L3.range.array.m[i] * cos(theta.radar.volume.L3.deg[k])
    }
  }
  
  # convert all tilt dBZ (log scale) to Ze (Z-effective, non-log scale)
  #   Ze is used to define the lowest height level classified by NSSL HCA category (the operational L3 product) category 'hail with rain')
  print(paste("  Convert L2 Z to Ze, for later calculations", sep=""))
  nexrad.L2.tilt1.Ze.new                  <- 10 ^ (nexrad.L2.tilt1.Z.new  / 10)
  nexrad.L2.tilt2.Ze.new                  <- 10 ^ (nexrad.L2.tilt2.Z.new  / 10)
  nexrad.L2.tilt3.Ze.new                  <- 10 ^ (nexrad.L2.tilt3.Z.new  / 10)
  nexrad.L2.tilt4.Ze.new                  <- 10 ^ (nexrad.L2.tilt4.Z.new  / 10)
  nexrad.L2.tilt5.Ze.new                  <- 10 ^ (nexrad.L2.tilt5.Z.new  / 10)
  nexrad.L2.tilt6.Ze.new                  <- 10 ^ (nexrad.L2.tilt6.Z.new  / 10)

  # load first 6 tilts of Ze.new into a 3D matrix
  nexrad.L2.tilt1to6.Ze.new               <- array(data = NA, dim = c(dim(nexrad.L2.tilt1.Z.new)[1], dim(nexrad.L2.tilt1.Z.new)[2], length(theta.radar.volume.L3.deg)))
  nexrad.L2.tilt1to6.Ze.new[,,1]          <- nexrad.L2.tilt1.Ze.new
  nexrad.L2.tilt1to6.Ze.new[,,2]          <- nexrad.L2.tilt2.Ze.new
  nexrad.L2.tilt1to6.Ze.new[,,3]          <- nexrad.L2.tilt3.Ze.new
  nexrad.L2.tilt1to6.Ze.new[,,4]          <- nexrad.L2.tilt4.Ze.new
  nexrad.L2.tilt1to6.Ze.new[,,5]          <- nexrad.L2.tilt5.Ze.new
  nexrad.L2.tilt1to6.Ze.new[,,6]          <- nexrad.L2.tilt6.Ze.new

  # load first 6 tilts of L3.HCA into a 3D matrix, with only 1 degree azimuth rez
  nexrad.L3.tilt1to6.HCA                  <- array(data = NA, dim = c(dim(nexrad.L3.tilt1.HCA)[1], dim(nexrad.L3.tilt1.HCA)[2], length(theta.radar.volume.L3.deg)))
  if (length(dim(nexrad.L3.tilt1.HCA))==3) {
    nexrad.L3.tilt1to6.HCA[,,1]             <- nexrad.L3.tilt1.HCA[,,2]
  } else {
    nexrad.L3.tilt1to6.HCA[,,1]             <- nexrad.L3.tilt1.HCA[,]
  }
  if (length(dim(nexrad.L3.tilt2.HCA))==3) {
    nexrad.L3.tilt1to6.HCA[,,2]             <- nexrad.L3.tilt2.HCA[,,2]
  } else {
    nexrad.L3.tilt1to6.HCA[,,2]             <- nexrad.L3.tilt2.HCA[,]
  }
  if (length(dim(nexrad.L3.tilt3.HCA))==3) {
    nexrad.L3.tilt1to6.HCA[,,3]             <- nexrad.L3.tilt3.HCA[,,2]
  } else {
    nexrad.L3.tilt1to6.HCA[,,3]             <- nexrad.L3.tilt3.HCA[,]
  }
  if (length(dim(nexrad.L3.tilt4.HCA))==3) {
    nexrad.L3.tilt1to6.HCA[,,4]             <- nexrad.L3.tilt4.HCA[,,2]
  } else {
    nexrad.L3.tilt1to6.HCA[,,4]             <- nexrad.L3.tilt4.HCA[,]
  }
  if (length(dim(nexrad.L3.tilt5.HCA))==3) {
    nexrad.L3.tilt1to6.HCA[,,5]             <- nexrad.L3.tilt5.HCA[,,2]
  } else {
    nexrad.L3.tilt1to6.HCA[,,5]             <- nexrad.L3.tilt5.HCA[,]
  }
  if (length(dim(nexrad.L3.tilt6.HCA))==3) {
    nexrad.L3.tilt1to6.HCA[,,6]             <- nexrad.L3.tilt6.HCA[,,2]
  } else {
    num.row.diff                            <- dim(nexrad.L3.tilt1.HCA)[1] - dim(nexrad.L3.tilt6.HCA)[1]
    nexrad.L3.tilt6.HCA                     <- rbind(nexrad.L3.tilt6.HCA[,], array(data = NA, dim = c(num.row.diff, dim(nexrad.L3.tilt1.HCA)[2])) )
    nexrad.L3.tilt1to6.HCA[,,6]             <- nexrad.L3.tilt6.HCA
  }
  
  # expand 3D matrix to 720 azimuths (0.5 degree rez)
  print(paste("  Replicate L3 HCA azim dim from 1.0 degree to 0.5 degree azim rez", sep=""))
  nexrad.L3.tilt1to6.HCA.new              <- array(data = NA, dim = c(dim(nexrad.L3.tilt1.HCA)[1], dim(nexrad.L2.tilt1.Z.new)[2], length(theta.radar.volume.L3.deg)))
  for (j in 1:dim(nexrad.L3.tilt1to6.HCA)[2]) {
    nexrad.L3.tilt1to6.HCA.new[,j*2,]     <- nexrad.L3.tilt1to6.HCA[,j,]
    nexrad.L3.tilt1to6.HCA.new[,(j*2)-1,] <- nexrad.L3.tilt1to6.HCA[,j,]
  }
  
  # find unique cell id values
  print(paste("  Find unique L3 SSTRUC cell IDs for volume #", a, sep=""))
  cell.id.unique                         <- unique(nexrad.L3.SSTRUC.alltimes.df$cell.id[ind.L3.SSTRUC.all.closest])
  #loop through and plot 
  for (t in 1:length(cell.id.unique)) {
    ind.cell.id <- which(nexrad.L3.SSTRUC.alltimes.df$cell.id == cell.id.unique[t])
    #print(t)
    #print(ind.cell.id)
    if (t == 1) {
      par(mar=c(4, 4, 3, 3))
      plot(nexrad.L3.SSTRUC.alltimes.df$lon.deg[ind.cell.id], nexrad.L3.SSTRUC.alltimes.df$lat.deg[ind.cell.id], col="cyan", xlim=c(min(nexrad.L3.SSTRUC.alltimes.df$lon.deg), max(nexrad.L3.SSTRUC.alltimes.df$lon.deg)), ylim=c(min(nexrad.L3.SSTRUC.alltimes.df$lat.deg), max(nexrad.L3.SSTRUC.alltimes.df$lat.deg)), type="b")
      grid()
    } else if (t == 8) {
      lines(nexrad.L3.SSTRUC.alltimes.df$lon.deg[ind.cell.id], nexrad.L3.SSTRUC.alltimes.df$lat.deg[ind.cell.id], col="red", type = "b")
      text(nexrad.L3.SSTRUC.alltimes.df$lon.deg[tail(ind.cell.id, n=1)]+0.2, nexrad.L3.SSTRUC.alltimes.df$lat.deg[tail(ind.cell.id, n=1)], cell.id.unique[t], col="red")
    } else {
      lines(nexrad.L3.SSTRUC.alltimes.df$lon.deg[ind.cell.id], nexrad.L3.SSTRUC.alltimes.df$lat.deg[ind.cell.id], col="cyan", type ="b")
    }
    text(radar.lon.deg, radar.lat.deg,     "*",    col="black")
    text(radar.lon.deg, radar.lat.deg+0.1, "KAMA", col="black")
    
  }    # end of for (t in 1:...)
  
  #plot(nexrad.L3.SSTRUC.alltimes.df$hgt.z.max.kft[8], type="l")

  #-----------------------------
  # equations (from Wallace and Friedrich, WAF, 2019)
  #-----------------------------
  print(paste("----------------------------------", sep=""))
  print(paste("Begin calculations for volume #", a, sep=""))
  print(paste("  Define values using temp profile", sep=""))
  
  # STARTING AT SFC, CALC Temp.dewpoint.depress AT EACH LEVEL UNTIL VALUE CROSSES ZERO
  print(paste("    Define dewpoint depress", sep=""))
  Temp.dewpoint.depress        <- matrix(0L, nrow = dim(balloon.sounding.df)[1], ncol = 1)
  Temp.dewpoint.depress[1]     <- NA
  index.dewpoint.depress       <- 0
  for (i in 2:dim(balloon.sounding.df)[1]) {
    # Dewpoint depression shortcut calculation used by forecasters, called '1/3 rule'
    Temp.dewpoint.depress[i] <- (balloon.sounding.df$temp.amb.C[i] - balloon.sounding.df$temp.dew.C[i]) / 3 
    if (Temp.dewpoint.depress[i] <= 0.2 & flag.first.time == 0) {
      index.dewpoint.depress <- i
      flag.first.time        <- 1
    }
  }    # end of for (i in 2:...)
  # IN MOST SOUNDINGS, DEWPOINT DEPRESS NEVER ACTUALLY GETS TO ZERO, FIND MIN VALUE
  if (index.dewpoint.depress == 0) {
    ind.dewpoint.depress     <- which.min(Temp.dewpoint.depress)
  }

  # DEFINE A NUMBER OF VALUES FROM THE SINGLE STATIC TEMPERATURE PROFILE
  #   calculate pressure(mean) from sfc to zero degree wet bulb calculation
  print(paste("    Define pres of 0 deg wetbulb and mean pres from sfc to 0 deg wetbulb", sep=""))
  Pres.sfc.mb                  <- balloon.sounding.df$Pres.hPa[1]
  Pres.0degwetbulb.mb          <- balloon.sounding.df$Pres.hPa[ind.dewpoint.depress]
  Pres.mean.sfcto0degwetbulb.mb<- (Pres.sfc.mb + Pres.0degwetbulb.mb) / 2

  #   define indices for required temperature values
  ind.first.subzero.temp.amb.C <- min(which(balloon.sounding.df$temp.amb.C < 0))
  ind.first.subm20.temp.amb.C  <- min(which(balloon.sounding.df$temp.amb.C < -20))
  #  H.ARL.sub0   is height above radar level of melting level, SINGLE VALUE FOR RADAR DOMAIN
  H.ARL.sub0                   <- balloon.sounding.df$hgt.m[ind.first.subzero.temp.amb.C] -  radar.elev.m  # [m]
  #  H.ARL.subm20 is height above radar level of the -20 deg C ambient temp, SINGLE VALUE FOR RADAR DOMAIN
  H.ARL.subm20                 <- balloon.sounding.df$hgt.m[ind.first.subm20.temp.amb.C]  -  radar.elev.m  # [m]
  print(paste("    Define H.ARL.sub0 and H.ARL.subm20", sep=""))
  
  # pressure correction, defined as mean pressure between sfc and height of 0 deg wet bulb temp from utilized sounding
  # CHECK THE PUB, ENSURE THIS IS A SCALE FACTOR [UNITLESS]
  print(paste("    Define pressure correction factor     [from mean sfc to 0 deg wetbulb pres]", sep=""))
  Pres.correction.factor                           <- (1000 * Pres.mean.sfcto0degwetbulb.mb ^ -1) ^ 0.545                                       # input a single scalar, in [mb]

  print(paste("  Define values using L2 and L3 data"))

  # IWC.ofhail.subtemp is calculated by only using Ze data below the height of the 0 deg isotherm (below H.ARL.sub0) and at the lowest height in a 'column' where RAIL/HAIL HCA exists
  #   from here and onward, all ranges confined to max range of L3 data = 1200 gates
  print(paste("    Define index where H.ARL is below FZLVL hgt [from sounding]", sep=""))
  ind.lt.H.ARL.sub0                                <- which(H.ARL.L2.m[1:1200,1,1:6] < H.ARL.sub0)
  print(paste("    Define index where L3 HCA is hail (val=9) or hail mixed with rain (val=10)", sep=""))
  ind.HCA.is.RAHA                                  <- which(nexrad.L3.tilt1to6.HCA.new[1:1200,,] == 90 | nexrad.L3.tilt1to6.HCA.new[1:1200,,] == 100)
  print(paste("    Define index where in both of these sets", sep=""))
  ind.in.both.sets                                 <- intersect(ind.lt.H.ARL.sub0, ind.HCA.is.RAHA)
  
  nexrad.L2.tilt1to6.Ze.bothsets                   <- array(data = 0L, dim=c(dim(nexrad.L3.tilt1.HCA)[1], dim(nexrad.L2.tilt1.Z.new)[2], length(theta.radar.volume.L3.deg)))
  nexrad.L2.tilt1to6.Ze.bothsets[ind.in.both.sets] <- nexrad.L2.tilt1to6.Ze.new[ind.in.both.sets]

  # Calculate Sev.Hail.Ind and MEHS.mm
  #   From here to Sev.Hail.Ing calculation is run for every 2D storm cell detected by a cell tracking algorithm (SCIT, Johnson et al., 1998)
  #   Z (which is dBZ) is the maximum refl value for each storm component, applied across vertical depth of storm (hgt.top.ft - hgt.base.ft)
  #   below is calculated using info from the 2D storm components for each cell being analyzed, need at least two components
  #   E calculated using the max refl for each storm component, then applied across the depth of the storm component
  nexrad.L3.tilt1.LAT720                           <- array(data = NA, dim = c(dim(nexrad.L3.tilt1.LAT)[1], 2*dim(nexrad.L3.tilt1.LAT)[2]))
  nexrad.L3.tilt1.LON720                           <- array(data = NA, dim = c(dim(nexrad.L3.tilt1.LON)[1], 2*dim(nexrad.L3.tilt1.LON)[2]))
  for (j in 1:dim(nexrad.L3.tilt1.LAT)[2]){
    nexrad.L3.tilt1.LAT720[,(j*2)]        <- nexrad.L3.tilt1.LAT[,j]
    nexrad.L3.tilt1.LAT720[,(j*2)-1]      <- nexrad.L3.tilt1.LAT[,j]
    nexrad.L3.tilt1.LON720[,(j*2)]        <- nexrad.L3.tilt1.LON[,j]
    nexrad.L3.tilt1.LON720[,(j*2)-1]      <- nexrad.L3.tilt1.LON[,j]
  }

  # print stuff to GUI
  print(paste("    Loop over all ", length(ind.L3.SSTRUC.all.closest), " L3 storm cells in volume #", a,      sep=""))
  print(paste("      Define W.subZ weighting         [from Zmax val compared to constants, 0-1 scale]",       sep=""))
  print(paste("      Calc kinetic energy of hail     [from Zmax val/weighting]",                              sep=""))
  print(paste("      Define W.subT.ofH weighting     [from Zmax hgt val compared to m20 hgt val, 0-1 scale]", sep=""))
  print(paste("      Calc delta H.subT               [from L3 SSTRUC hgt top-base]",                          sep=""))
  print(paste("      Calc Sev Hail Ind               [from k-energy/weighting/deltaH.subT]",                  sep=""))
  print(paste("      Calc MEHS                       [from Sev Hail Index val]",                              sep=""))
  print(paste("      Calc velocity of hail           [from MEHS/pres correction factor]",                     sep=""))
  print(paste("      Calc sfc accum of hail for 1vol [from constants/IWC-of-hail/velocity-of-hail]",          sep=""))
  
  # initialize arrays to be filled in following loop
  Energy.kinetic.ofHail                            <- array(data = 0L, dim = c(dim(nexrad.L3.SSTRUC.alltimes.df)[1]))
  X.rad                                            <- array(data = 0L, dim = c(dim(nexrad.L3.SSTRUC.alltimes.df)[1]))
  Y.rad                                            <- array(data = 0L, dim = c(dim(nexrad.L3.SSTRUC.alltimes.df)[1]))
  d                                                <- array(data = 0L, dim = c(dim(nexrad.L3.SSTRUC.alltimes.df)[1]))
  c                                                <- array(data = 0L, dim = c(dim(nexrad.L3.SSTRUC.alltimes.df)[1]))
  bearing.radartoSSTRUC.deg                        <- array(data = NA, dim = c(dim(nexrad.L3.SSTRUC.alltimes.df)[1]))
  dist.radartoSSTRUC.km                            <- array(data = NA, dim = c(dim(nexrad.L3.SSTRUC.alltimes.df)[1]))
  W.subZ                                           <- array(data = 0L, dim = c(dim(nexrad.L3.SSTRUC.alltimes.df)[1]))
  W.subT.ofH                                       <- array(data = 0L, dim = c(dim(nexrad.L3.SSTRUC.alltimes.df)[1]))
  delta.H.subT                                     <- array(data = 0L, dim = c(dim(nexrad.L3.SSTRUC.alltimes.df)[1]))
  Sev.Hail.Ind.tot                                 <- array(data = 0L, dim = c(dim(nexrad.L3.SSTRUC.alltimes.df)[1]))
  MEHS.cm                                          <- array(data = 0L, dim = c(dim(nexrad.L3.SSTRUC.alltimes.df)[1]))
  Vel.t.of.diam                                    <- array(data = 0L, dim = c(dim(nexrad.L3.SSTRUC.alltimes.df)[1]))
  ind.L3.tilt1.LATLON.closesttocell                <- array(data = 0L, dim = c(dim(nexrad.L3.SSTRUC.alltimes.df)[1]))
  ind.min.dist.stormcell.to.each.point             <- array(data = 0L, dim = c(dim(nexrad.L3.SSTRUC.alltimes.df)[1]))
  IWC.ofHail.subtemp.gm3                           <- array(data = 0L, dim = c(dim(nexrad.L3.SSTRUC.alltimes.df)[1]))
  sfc.hail.1volaccum.cm                            <- array(data = 0L, dim = c(dim(nexrad.L3.SSTRUC.alltimes.df)[1]))
  # loop over all storm cells identified by the L3 SSTRUC product at this volume time
  for (m in 1:length(ind.L3.SSTRUC.all.closest)) {
  #for (m in 1:dim(nexrad.L3.SSTRUC.alltimes.df)[1]) {
  
    # calculate bearing (degrees) between two lat/lon points
    #β = atan2(X,Y) is bearing from X to Y
    #X = cos θb * sin ∆L
    #Y = cos θa * sin θb – sin θa * cos θb * cos ∆L
    X.rad[m]                          <- cos(NISTdegTOradian(nexrad.L3.SSTRUC.alltimes.df$lat.deg[m])) * sin(abs(NISTdegTOradian(radar.lon.deg)-NISTdegTOradian(nexrad.L3.SSTRUC.alltimes.df$lon.deg[m])))
    Y.rad[m]                          <- cos(NISTdegTOradian(radar.lat.deg)) * sin(NISTdegTOradian(nexrad.L3.SSTRUC.alltimes.df$lat.deg[m])) - sin(NISTdegTOradian(radar.lat.deg)) * cos(NISTdegTOradian(nexrad.L3.SSTRUC.alltimes.df$lat[m])) * cos(abs(NISTdegTOradian(radar.lon.deg)-NISTdegTOradian(nexrad.L3.SSTRUC.alltimes.df$lon[m])))
    bearing.radartoSSTRUC.deg[m]      <- 360 - NISTradianTOdeg(atan2(X.rad[m], Y.rad[m]))
    # find ind of closest bearing to a L2 azim
    ind.closest.azim                  <- which.min(abs(nexrad.L2.tilt1.AZIM-bearing.radartoSSTRUC.deg[m]))
    #calculate great circle distance along earth's surface between two lat/lon pts
    #a = sin²(Δφ/2) + cos φ1 ⋅ cos φ2 ⋅ sin²(Δλ/2)
    #c = 2 ⋅ atan2( √a, √(1−a) )
    #d = R ⋅ c
    d[m]                              <- sin(abs(NISTdegTOradian(radar.lat.deg) - NISTdegTOradian(nexrad.L3.SSTRUC.alltimes.df$lat.deg[m])) / 2) * sin(abs(NISTdegTOradian(radar.lat.deg)-NISTdegTOradian(nexrad.L3.SSTRUC.alltimes.df$lat.deg[m])) / 2) + (cos(NISTdegTOradian(radar.lat.deg)) * cos(NISTdegTOradian(nexrad.L3.SSTRUC.alltimes.df$lat.deg[m])) * sin(abs(NISTdegTOradian(radar.lon.deg) - NISTdegTOradian(nexrad.L3.SSTRUC.alltimes.df$lon.deg[m])) / 2) * sin(abs(NISTdegTOradian(radar.lon.deg)-NISTdegTOradian(nexrad.L3.SSTRUC.alltimes.df$lon.deg[m])) / 2) )
    c[m]                              <- 2 * atan2(sqrt(d[m]), sqrt(1-d[m]))
    dist.radartoSSTRUC.km[m]          <- radius.earth.km * c[m]
    #
    ind.closest.horiz.dist            <- round((which.min(abs(dist.radartoSSTRUC.km[m] - Dist.from.radar.m[,1,1]/1000))) / 2)
    
    # calc IWC of Hail
    #print(paste("    Calc IWC of hail 3D matrix where in both of these sets"))
    IWC.ofHail.subtemp.gm3[m]         <- (4.4 ^ 10e-5) * (mean(nexrad.L2.tilt1.Ze.new[seq(from = ind.closest.horiz.dist-5, to = ind.closest.horiz.dist+5, by = 1), ind.closest.azim], na.rm=T) ^ 0.71)  # [g/m3], as per Kalina et al, 2016, input 3D matrix
    #print(IWC.ofHail.subtemp.gm3[m])
    
    # define weight W.subZ for L3 SSTRUC storm cell number m
    if (nexrad.L3.SSTRUC.alltimes.df$z.max.dbz[ind.L3.SSTRUC.all.closest[m]] <= Z.subL | is.na(nexrad.L3.SSTRUC.alltimes.df$z.max.dbz[ind.L3.SSTRUC.all.closest[m]])) {
      W.subZ[m] <- 0
    } else if (nexrad.L3.SSTRUC.alltimes.df$z.max.dbz[ind.L3.SSTRUC.all.closest[m]] >= Z.subU) {
      W.subZ[m] <- 1
    } else if (nexrad.L3.SSTRUC.alltimes.df$z.max.dbz[ind.L3.SSTRUC.all.closest[m]] > Z.subL & nexrad.L3.SSTRUC.alltimes.df$z.max.dbz[ind.L3.SSTRUC.all.closest[m]] < Z.subU) {
      W.subZ[m] <- (nexrad.L3.SSTRUC.alltimes.df$z.max.dbz[ind.L3.SSTRUC.all.closest[m]] - Z.subL)/(Z.subU - Z.subL)
    }  # end of if (z.max.dbz[m] <= Z.subL)

    Energy.kinetic.ofHail[m]                <- (5 * 10 ^ -6) * (10 ^ (0.084 * nexrad.L3.SSTRUC.alltimes.df$z.max.dbz[ind.L3.SSTRUC.all.closest[m]])) * W.subZ[m]                          # [J m-2 s-1], both inputs 3D matrices
    
    # LOOP OVER NUMBER OF LEVELS OF STORM, WHERE EACH LEVEL REFERS TO HOW MANY TILTS PASS THROUGH STORM FROM HSUB0 TO HSUBTOP
    #  ENERGY.KINETIC.OFHAIL IS NOW A CONSTANT FOR EACH STORM COMPONENT
    #  W.SUBT.OFH AND DELTA.H.SUBT ARE CALCULATED AT EACH LEVEL
    #############################
    # CALCULATE NUM TILT LEVELS PASS THROUGH EACH STORM CENTROID
    #############################
    #num.tilts.through.storm[m]    <- 4
    #for (n in 1:num.tilts.through.storm[m]) {
    if (nexrad.L3.SSTRUC.alltimes.df$hgt.z.max.kft[ind.L3.SSTRUC.all.closest[m]]*1000 <= H.ARL.sub0) {
      W.subT.ofH[m] <- 0
    } else if (nexrad.L3.SSTRUC.alltimes.df$hgt.z.max.kft[ind.L3.SSTRUC.all.closest[m]]*1000 >= H.ARL.subm20) {
      W.subT.ofH[m] <- 1
    } else if (nexrad.L3.SSTRUC.alltimes.df$hgt.z.max.kft[ind.L3.SSTRUC.all.closest[m]]*1000 > H.ARL.sub0 & nexrad.L3.SSTRUC.alltimes.df$hgt.z.max.kft[ind.L3.SSTRUC.all.closest[m]]*1000 < H.ARL.subm20) {
      W.subT.ofH[m] <- (nexrad.L3.SSTRUC.alltimes.df$hgt.z.max.kft[ind.L3.SSTRUC.all.closest[m]]*1000 - H.ARL.sub0) / (H.ARL.subm20 - H.ARL.sub0)
    }  # end of if (hgt.z.max.kft[m] <= H.ARL.sub0)
    
    delta.H.subT[m]                         <- nexrad.L3.SSTRUC.alltimes.df$hgt.top.kft[ind.L3.SSTRUC.all.closest[m]] - nexrad.L3.SSTRUC.alltimes.df$hgt.base.kft[ind.L3.SSTRUC.all.closest[m]]
    Sev.Hail.Ind.tot[m]                     <- 0.1 * (W.subT.ofH[m] * Energy.kinetic.ofHail[m] * delta.H.subT[m]) # [J m-1 s-1], from Witt et al., 1998, performs well S. Plains, weaker for further east
    #Sev.Hail.Ind.tot[m]                  <- Sev.Hail.Ind.tot + Sev.Hail.Ind.tilt
    #}  # end of for (n in 1:num.tilts.through.storm)

    # max estimated hail size, a diameter (D) [cm]
    MEHS.cm[m]                              <- 2.54 * (Sev.Hail.Ind.tot[m] ^ 0.5)                                                                             # in [mm]
  
    # fall velocity relationship as a function of diam[cm], denoted as HW14Ri (Heymsfield and Wright, 2014, from Rimed crystals)
    Vel.t.of.diam[m]                        <- 488 * (MEHS.cm[m] ^ 0.84) * Pres.correction.factor * mPERcm    # input MEHS.mm [mm] is a 3D matrix, other a unitless scalar
  
    # find lat/lon from L3 SSTRUC that are spatially closest to cell m
    L3.tilt1.LAT720.df                      <- data.frame(as.vector(nexrad.L3.tilt1.LAT720))
    L3.tilt1.LON720.df                      <- data.frame(as.vector(nexrad.L3.tilt1.LON720))
    L3.tilt1.LATLON720.df                   <- cbind(L3.tilt1.LON720.df, L3.tilt1.LAT720.df)
    L3.SSTRUC.LATLON.df                     <- data.frame(c(nexrad.L3.SSTRUC.alltimes.df$lat.deg[m], nexrad.L3.SSTRUC.alltimes.df$lon.deg[m]))
    
    ## THIS TOO TAKES WAY TOO LONG
    #dist.stormcell.to.each.point            <- array(data = 0L, dim = c(length(L3.tilt1.LAT720.df)))
    #for (a in 1:dim(L3.tilt1.LAT720.df)[1]) {
    #  #if (a/1000 print(a)
    #  dist.stormcell.to.each.point[a] <- distm(c(L3.SSTRUC.LATLON.df$c.nexrad.L3.SSTRUC.alltimes.df.lat.deg.m...nexrad.L3.SSTRUC.alltimes.df.lon.deg.m..[2], L3.SSTRUC.LATLON.df$c.nexrad.L3.SSTRUC.alltimes.df.lat.deg.m...nexrad.L3.SSTRUC.alltimes.df.lon.deg.m..[1]), c(L3.tilt1.LON720.df$as.vector.nexrad.L3.tilt1.LON720.[a], L3.tilt1.LAT720.df$as.vector.nexrad.L3.tilt1.LAT720.[a]), fun = distVincentyEllipsoid)
    #}
    #ind.min.dist.stormcell.to.each.point[m] <- which.min(dist.stormcell.to.each.point)
    #
    #IWC.ofHail.subtemp.cellarray.gm3[m]     <- IWC.ofHail.subtemp.tilt1.gm3[ind.min.dist.stormcell.to.each.point[m]] #
    
    
    # Calculate accumulation of hail depth [cm] at the surface
    #   set up loop over each radar volume (time) to sum the accumulation [cm] at each radar pixel [range, azimuth]
    #   accum.sfc.hail is the depth accum within the time.delta that represents the given radar volume ( in this case, 4.1 minutes)
    
    # [3D MATRIX]                           <-     [scalar]         [scalar]          [scalar]     [scalar]    [3D matrix]                          [3D matrix]               [scalar]   [scalar]        [scalar]
    # [UNITS=cm]                            <-     [NA]             [1/NA]            [cm3/g]      [m3/cm3]    [g/m3]                               [m/s]                     [cm/m]     [min]           [s/min]
    #sfc.hail.1volaccum.cm[radar.volume]    <- (1 / epsilon) * (1 / packing.density * rho.hail) * (m3PERcm3) * (IWC.ofHail.subtemp.gm3[1:1200,,] * (Vel.t.of.diam[1:1200,,] * cmPERm) * (time.delta.m * secPERmin))
    sfc.hail.1volaccum.cm[m]                <- (1 / epsilon) * (1 / packing.density * rho.hail) * (m3PERcm3) * (IWC.ofHail.subtemp.gm3[m] * (Vel.t.of.diam[m] * cmPERm) * (time.delta.m * secPERmin))
    
    # accumulate sfc.hail.1volaccum.cm values around each m storm cell 
    #ind.seq.horiz.dist                                               <- seq(from = ind.closest.horiz.dist-5, to = ind.closest.horiz.dist+5, by = 1)
    #ind.seq.azim                                                     <- seq(from = ind.closest.azim-4, to = ind.closest.azim+4, by = 1)
    #sfc.hail.1volaccum.cm[ind.seq.horiz.dist, round(ind.seq.azim/2)] <- sfc.hail.1volaccum.cm[m]
    #sfc.hail.totaccum.cm                                             <- sfc.hail.totaccum.cm + sfc.hail.1volaccum.cm
    
    # print some values to GUI
    print(paste("        for storm cell #", m, " of ", length(ind.L3.SSTRUC.all.closest), sep=""))
    print(paste("          ", nexrad.L3.SSTRUC.alltimes.df$cell.id[ind.L3.SSTRUC.all.closest[m]], " lat/lon: ", nexrad.L3.SSTRUC.alltimes.df$lat.deg[ind.L3.SSTRUC.all.closest[m]], ", ", nexrad.L3.SSTRUC.alltimes.df$lon.deg[ind.L3.SSTRUC.all.closest[m]], sep=""))
    print(paste("             top: ", nexrad.L3.SSTRUC.alltimes.df$hgt.top.kft[ind.L3.SSTRUC.all.closest[m]], " base: ", nexrad.L3.SSTRUC.alltimes.df$hgt.base.kft[ind.L3.SSTRUC.all.closest[m]], " Zmax: ", nexrad.L3.SSTRUC.alltimes.df$z.max.dbz[ind.L3.SSTRUC.all.closest[m]], sep=""))
    print(paste("             bearing:                                  = ", bearing.radartoSSTRUC.deg[m]),                                                        sep="")
    #print(paste("             multiply all constants                    = ", (1/epsilon)*(1/packing.density*rho.hail)*(m3PERcm3)*cmPERm*(time.delta.m*secPERmin)), sep="")
    print(paste("             IWC.ofHail.subtemp.gm3[", m, "]           = ", IWC.ofHail.subtemp.gm3[m]),                                                           sep="")
    print(paste("             Vel.t.of.diam[", m, "]                    = ", Vel.t.of.diam[m]),                                                                    sep="")                                                
    print(paste("             sfc.hail.1volaccum.cm[", m, "]            = ", sfc.hail.1volaccum.cm[m]),                                                            sep="")
    
    #sfc.hail.totaccum.cm     <- sfc.hail.totaccum.cm + sfc.hail.1volaccum.cm                                                      # matrix addition

  }    # end of for (m in 1:dim(nexrad.L3.SSTRUC.0530Z.df)[1]) 
    
  # APPLY SFC HAIL ACCUM FOR ALL m STORM CELLS
  print(paste("  Build data frame of all calculated fields for volume #", a, sep=""))
  
  if (exists("nexrad.L3.SSTRUC.1vol.df")) {
    rm(nexrad.L3.SSTRUC.1vol.df)
  }
  
  # Add MEHS.cm to SSTRUC df for plotting
  MEHS.cm.df                            <- data.frame(MEHS.cm)
  MEHS.cm.1vol.df                       <- data.frame(MEHS.cm[1:7])
  nexrad.L3.SSTRUC.1vol.df              <- cbind(nexrad.L3.SSTRUC.alltimes.df[1:7,], MEHS.cm.1vol.df)
  nexrad.L3.SSTRUC.alltimes.df          <- cbind(nexrad.L3.SSTRUC.alltimes.df, MEHS.cm.df)

  # Add Vel.t.of.diam to SSTRUC df for plotting
  Vel.t.of.diam.df                      <- data.frame(Vel.t.of.diam)
  Vel.t.of.diam.1vol.df                 <- data.frame(Vel.t.of.diam[1:7])
  nexrad.L3.SSTRUC.alltimes.df          <- cbind(nexrad.L3.SSTRUC.alltimes.df, Vel.t.of.diam.df)
  nexrad.L3.SSTRUC.1vol.df              <- cbind(nexrad.L3.SSTRUC.1vol.df, Vel.t.of.diam.1vol.df)

  # Add IWC.ofHAIL.subtemp.cellarray.gm3 to SSTRUC df for plotting
  IWC.ofHail.subtemp.gm3.df             <- data.frame(IWC.ofHail.subtemp.gm3)
  IWC.ofHail.subtemp.gm3.1vol.df        <- data.frame(IWC.ofHail.subtemp.gm3[1:7])
  nexrad.L3.SSTRUC.alltimes.df          <- cbind(nexrad.L3.SSTRUC.alltimes.df, IWC.ofHail.subtemp.gm3.df)
  nexrad.L3.SSTRUC.1vol.df              <- cbind(nexrad.L3.SSTRUC.1vol.df, IWC.ofHail.subtemp.gm3.1vol.df)

  # Add sfc.hail.1volaccum.cm to SSTRUC df for plotting
  sfc.hail.1volaccum.cm.df              <- data.frame(sfc.hail.1volaccum.cm)
  sfc.hail.1volaccum.cm.1vol.df         <- data.frame(sfc.hail.1volaccum.cm[1:7])
  nexrad.L3.SSTRUC.alltimes.df          <- cbind(nexrad.L3.SSTRUC.alltimes.df, sfc.hail.1volaccum.cm.df)
  nexrad.L3.SSTRUC.1vol.df              <- cbind(nexrad.L3.SSTRUC.1vol.df, sfc.hail.1volaccum.cm.1vol.df)
  
  head(nexrad.L3.SSTRUC.alltimes.df)
  head(nexrad.L3.SSTRUC.1vol.df)
  
  #--------------------------------
  # Plotting
  #--------------------------------
  
  if (flag.make.plots == 1) {
  # area map for radar domain
  KAMA.map                            <- get_map(location = c(radar.lon.deg, radar.lat.deg), maptype = "terrain", source = "google", zoom = 8)

  print(paste("Begin plotting", sep=""))
  # plot Diam-max of hail
  ggmap(KAMA.map) + 
    geom_point(data=nexrad.L3.SSTRUC.1vol.df, aes(x=lon.deg, y=lat.deg),                 color="red",       size=1) +
    geom_point(data=nexrad.L3.SSTRUC.1vol.df, aes(x=lon.deg, y=lat.deg),                 color="darkgreen", size=1) +
    geom_point(data=nexrad.L3.SSTRUC.1vol.df, aes(x=lon.deg, y=lat.deg,                  color=MEHS.cm.1.7.),    size=5) +
    geom_text(data=nexrad.L3.SSTRUC.1vol.df,  aes(x=lon.deg, y=lat.deg,  label=cell.id),                    size=4, vjust = 0, nudge_y = 0.1) +
    geom_point(data=NEXRAD.site.df,           aes(x=lon,     y=lat),                                        size=4, shape=8) +
    geom_text(data=NEXRAD.site.df,            aes(x=lon,     y=lat,      label=ICAO),                       size=4, vjust = 0, nudge_y = 0.1) +
    labs(x = 'Lon', y = 'Lat') + 
    ggtitle(paste(substr(nexrad.L2.listfiles[1+a*b+(b-(b-0))], 12, 19), " @ ", substr(nexrad.L2.listfiles[1+a*b+(b-(b-0))], 21, 24), " GMT: Max diam. of hail [cm]", sep=""))

  # plot Vel
  ggmap(KAMA.map) + 
    geom_point(data=nexrad.L3.SSTRUC.1vol.df, aes(x=lon.deg, y=lat.deg),                 color="red",          size=1) +
    geom_point(data=nexrad.L3.SSTRUC.1vol.df, aes(x=lon.deg, y=lat.deg),                 color="darkgreen",    size=1) +
    geom_point(data=nexrad.L3.SSTRUC.1vol.df, aes(x=lon.deg, y=lat.deg,                  color=Vel.t.of.diam.1.7.), size=5) +
    geom_text(data=nexrad.L3.SSTRUC.1vol.df,  aes(x=lon.deg, y=lat.deg,  label=cell.id),                       size=4, vjust = 0, nudge_y = 0.1) +
    geom_point(data=NEXRAD.site.df,               aes(x=lon,     y=lat),                                           size=4, shape=8) +
    geom_text(data=NEXRAD.site.df,                aes(x=lon,     y=lat,      label=ICAO),                          size=4, vjust = 0, nudge_y = 0.1) +
    labs(x = 'Lon', y = 'Lat') + 
    ggtitle(paste(substr(nexrad.L2.listfiles[1+a*b+(b-(b-0))], 12, 19), " @ ", substr(nexrad.L2.listfiles[1+a*b+(b-(b-0))], 21, 24), " GMT: Velocity of hail [ms-1]", sep=""))

  # plot IWC.ofHail
  ggmap(KAMA.map) + 
    geom_point(data=nexrad.L3.SSTRUC.1vol.df, aes(x=lon.deg, y=lat.deg),                 color="red",          size=1) +
    geom_point(data=nexrad.L3.SSTRUC.1vol.df, aes(x=lon.deg, y=lat.deg),                 color="darkgreen",    size=1) +
    geom_point(data=nexrad.L3.SSTRUC.1vol.df, aes(x=lon.deg, y=lat.deg,                  color=IWC.ofHail.subtemp.gm3.1.7.), size=5) +
    geom_text(data=nexrad.L3.SSTRUC.1vol.df,  aes(x=lon.deg, y=lat.deg,  label=cell.id),                       size=4, vjust = 0, nudge_y = 0.1) +
    geom_point(data=NEXRAD.site.df,               aes(x=lon,     y=lat),                                           size=4, shape=8) +
    geom_text(data=NEXRAD.site.df,                aes(x=lon,     y=lat,      label=ICAO),                          size=4, vjust = 0, nudge_y = 0.1) +
    labs(x = 'Lon', y = 'Lat') + 
    ggtitle(paste(substr(nexrad.L2.listfiles[1+a*b+(b-(b-0))], 12, 19), " @ ", substr(nexrad.L2.listfiles[1+a*b+(b-(b-0))], 21, 24), " GMT: IWC of hail [gm-3]", sep=""))

  # plot sfc.hail.1volaccum.cm
  ggmap(KAMA.map) + 
    geom_point(data=nexrad.L3.SSTRUC.1vol.df, aes(x=lon.deg, y=lat.deg),                 color="red",                  size=1) +
    geom_point(data=nexrad.L3.SSTRUC.1vol.df, aes(x=lon.deg, y=lat.deg),                 color="darkgreen",            size=1) +
    geom_point(data=nexrad.L3.SSTRUC.1vol.df, aes(x=lon.deg, y=lat.deg,                  color=sfc.hail.1volaccum.cm.1.7.), size=5) +
    #geom_text(data=nexrad.L3.SSTRUC.alltimes.df,  aes(x=lon.deg, y=lat.deg,  label=cell.id),                               size=4, vjust = 0, nudge_y = 0.1) +
    geom_point(data=NEXRAD.site.df,               aes(x=lon,     y=lat),                                                   size=4, shape=8) +
    geom_text(data=NEXRAD.site.df,                aes(x=lon,     y=lat,      label=ICAO),                                  size=4, vjust = 0, nudge_y = 0.1) +
    labs(x = 'Lon', y = 'Lat') + 
    ggtitle(paste(substr(nexrad.L2.listfiles[1+a*b+(b-(b-0))], 12, 19), " @ ", substr(nexrad.L2.listfiles[1+a*b+(b-(b-0))], 21, 24), " GMT: Sfc hail accum. [cm]", sep=""))
  
  # plot sfc.hail.totaccum.cm for first test volume
  #levelplot(sfc.hail.totaccum.cm)
  
  ## test plotting
  #plot(array.num.azimuths.L2.Zge0[,1])
  #lines(array.num.azimuths.L2.Zge0[,2])
  #lines(array.num.azimuths.L2.Zge0[,8])
  
  # more test plotting
  #par(mfrow=c(1,1))
  #par(mar=c(4, 4, 3, 3))
  #image(nexrad.L2.tilt4.Z)
  #image(nexrad.L2.tilt4.Z.new)
  #nexrad.L2.tilt4.Z[,30]

  # plots for debugging
  #   all HCA values plotted (1-10 or 1-14)
  #levelplot(nexrad.L3.tilt1.HCA[,,2])
  #levelplot(nexrad.L3.tilt2.HCA[,,2])
  #levelplot(nexrad.L3.tilt3.HCA[,,2])
  #levelplot(nexrad.L3.tilt4.HCA[,])
  #levelplot(nexrad.L3.tilt5.HCA[,,2])
  #levelplot(nexrad.L3.tilt6.HCA[,,2])
  ##   only HCA = 9 or 10 plotted, all else set to HCA = 0
  #ind.HCA.not9or10                                <- which(nexrad.L3.tilt1.HCA[,,2] < 90 | nexrad.L3.tilt1.HCA[,,2] > 100)
  #nexrad.L3.tilt1.HCA.9or10only                   <- nexrad.L3.tilt1.HCA[,,2]
  #nexrad.L3.tilt1.HCA.9or10only[ind.HCA.not9or10] <- 0
  #levelplot(nexrad.L3.tilt1.HCA.9or10only/10, xlab="x-range [ ]", ylab="y-range [ ]", xaxt="No")
  
  } else {
    print("Plotting flag is turned off")
  }  # end of if (flag.make.plots) ...
  
}    # end of for (a in volume.num) ...

  ###########################################################
  # All this stuff below was replaced by using L3 SSTRUC inputs
  ###########################################################
  # this code replaced by use of SSTRUC storm cell data, which have Zmax and Zmax hgt values
  ## find the max H.ARL.L2.m that a Z value exists within the L2 radar volume
  ##  H.ARL.subtop is height above radar level of the top of the storm, SINGLE VALUE FOR RADAR DOMAIN
  ##  NOTE: ECHO TOP HEIGHT IN L3 data defined as Z VALUE GREATER OF EQUAL TO 18 dBZ
  ##H.ARL.subtop 
  #array.num.azimuths.L2.Zge0   <- matrix(0L, nrow = dim(nexrad.L2.tilt1.Z.new)[1], ncol = length(theta.radar.volume.L2.deg))
  #for (kk in 1:dim(nexrad.L2.tilt1.Z)[1]) {
  #  ind.L2.tilt1                      <- which(nexrad.L2.tilt1.Z.new[kk,] >= Z.value.for.echo.top)
  #  ind.L2.tilt1
  #  array.num.azimuths.L2.Zge0[kk,1]  <- length(ind.L2.tilt1)
  #}
  #for (kk in 1:dim(nexrad.L2.tilt2.Z)[1]) {
  #  ind.L2.tilt2                      <- which(nexrad.L2.tilt2.Z.new[kk,] >= Z.value.for.echo.top)
  #  ind.L2.tilt2
  #  array.num.azimuths.L2.Zge0[kk,2]  <- length(ind.L2.tilt2)
  #}
  #for (kk in 1:dim(nexrad.L2.tilt3.Z)[1]) {
  #  ind.L2.tilt3                      <- which(nexrad.L2.tilt3.Z.new[kk,] >= Z.value.for.echo.top)
  #  array.num.azimuths.L2.Zge0[kk,3]  <- length(ind.L2.tilt3)
  #}
  #for (kk in 1:dim(nexrad.L2.tilt4.Z)[1]) {
  #  ind.L2.tilt4                      <- which(nexrad.L2.tilt4.Z.new[kk,] >= Z.value.for.echo.top)
  #  array.num.azimuths.L2.Zge0[kk,4]  <- length(ind.L2.tilt4)
  #}
  #for (kk in 1:dim(nexrad.L2.tilt5.Z)[1]) {
  #  ind.L2.tilt5                      <- which(nexrad.L2.tilt5.Z.new[kk,] >= Z.value.for.echo.top)
  #  array.num.azimuths.L2.Zge0[kk,5]  <- length(ind.L2.tilt5)
  #}
  #for (kk in 1:dim(nexrad.L2.tilt6.Z)[1]) {
  #  ind.L2.tilt6                      <- which(nexrad.L2.tilt6.Z.new[kk,] >= Z.value.for.echo.top)
  #  array.num.azimuths.L2.Zge0[kk,6]  <- length(ind.L2.tilt6)
  #}
  #for (kk in 1:dim(nexrad.L2.tilt7.Z)[1]) {
  #  ind.L2.tilt7                      <- which(nexrad.L2.tilt7.Z.new[kk,] >= Z.value.for.echo.top)
  #  array.num.azimuths.L2.Zge0[kk,7]  <- length(ind.L2.tilt7)
  #}
  #for (kk in 1:dim(nexrad.L2.tilt8.Z)[1]) {  
  #  ind.L2.tilt8                      <- which(nexrad.L2.tilt8.Z.new[kk,] >= Z.value.for.echo.top)
  #  array.num.azimuths.L2.Zge0[kk,8]  <- length(ind.L2.tilt8)
  #}
  #for (kk in 1:dim(nexrad.L2.tilt9.Z)[1]) {
  #  ind.L2.tilt9                      <- which(nexrad.L2.tilt9.Z.new[kk,] >= Z.value.for.echo.top)
  #  array.num.azimuths.L2.Zge0[kk,9]  <- length(ind.L2.tilt9)
  #}
  #for (kk in 1:dim(nexrad.L2.tilt10.Z)[1]) {
  #  ind.L2.tilt10                     <- which(nexrad.L2.tilt10.Z.new[kk,] >= Z.value.for.echo.top)
  #  array.num.azimuths.L2.Zge0[kk,10] <- length(ind.L2.tilt10)
  #}
  #for (kk in 1:dim(nexrad.L2.tilt11.Z)[1]) {
  #  ind.L2.tilt11                     <- which(nexrad.L2.tilt11.Z.new[kk,] >= Z.value.for.echo.top)
  #  array.num.azimuths.L2.Zge0[kk,11] <- length(ind.L2.tilt11)
  #}
  #for (kk in 1:dim(nexrad.L2.tilt12.Z)[1]) {
  #  ind.L2.tilt12                     <- which(nexrad.L2.tilt12.Z.new[kk,] >= Z.value.for.echo.top)
  #  array.num.azimuths.L2.Zge0[kk,12] <- length(ind.L2.tilt12)
  #}
  #for (kk in 1:dim(nexrad.L2.tilt13.Z)[1]) {
  #  ind.L2.tilt13                     <- which(nexrad.L2.tilt13.Z.new[kk,] >= Z.value.for.echo.top)
  #  array.num.azimuths.L2.Zge0[kk,13] <- length(ind.L2.tilt13)
  #}
  #for (kk in 1:dim(nexrad.L2.tilt14.Z)[1]) {
  #  ind.L2.tilt14                     <- which(nexrad.L2.tilt14.Z.new[kk,] >= Z.value.for.echo.top)
  #  array.num.azimuths.L2.Zge0[kk,14] <- length(ind.L2.tilt14)
  #}     # end of for (kk)
  ## find an index of the maximum range that has more than 1 value greater than 18 dBZ for each tilt
  #ind.L2.tilt1.hasvalue  <- max(which(array.num.azimuths.L2.Zge0[, 1]  > 1))
  #ind.L2.tilt2.hasvalue  <- max(which(array.num.azimuths.L2.Zge0[, 2]  > 1))
  #ind.L2.tilt3.hasvalue  <- max(which(array.num.azimuths.L2.Zge0[, 3]  > 1))
  #ind.L2.tilt4.hasvalue  <- max(which(array.num.azimuths.L2.Zge0[, 4]  > 1))
  #ind.L2.tilt5.hasvalue  <- max(which(array.num.azimuths.L2.Zge0[, 5]  > 1))
  #ind.L2.tilt6.hasvalue  <- max(which(array.num.azimuths.L2.Zge0[, 6]  > 1))
  #ind.L2.tilt7.hasvalue  <- max(which(array.num.azimuths.L2.Zge0[, 7]  > 1))
  #ind.L2.tilt8.hasvalue  <- max(which(array.num.azimuths.L2.Zge0[, 8]  > 1))
  #ind.L2.tilt9.hasvalue  <- max(which(array.num.azimuths.L2.Zge0[, 9]  > 1))
  #ind.L2.tilt10.hasvalue <- max(which(array.num.azimuths.L2.Zge0[, 10] > 1))
  #ind.L2.tilt11.hasvalue <- max(which(array.num.azimuths.L2.Zge0[, 11] > 1))
  #ind.L2.tilt12.hasvalue <- max(which(array.num.azimuths.L2.Zge0[, 12] > 1))
  #ind.L2.tilt13.hasvalue <- max(which(array.num.azimuths.L2.Zge0[, 13] > 1))
  #ind.L2.tilt14.hasvalue <- max(which(array.num.azimuths.L2.Zge0[, 14] > 1))
  ## create an array of the highest echo heights observed within each 2 tilt
  #H.ARL.max.fullvolume.L2.array.m <- c(H.ARL.L2.m[ind.L2.tilt1.hasvalue, , 1], H.ARL.L2.m[ind.L2.tilt2.hasvalue, , 2], H.ARL.L2.m[ind.L2.tilt3.hasvalue, , 3], H.ARL.L2.m[ind.L2.tilt4.hasvalue, , 4], H.ARL.L2.m[ind.L2.tilt5.hasvalue, , 5], H.ARL.L2.m[ind.L2.tilt6.hasvalue, , 6], H.ARL.L2.m[ind.L2.tilt7.hasvalue, , 7], H.ARL.L2.m[ind.L2.tilt8.hasvalue, , 8], H.ARL.L2.m[ind.L2.tilt9.hasvalue, , 9], H.ARL.L2.m[ind.L2.tilt10.hasvalue, , 10], H.ARL.L2.m[ind.L2.tilt11.hasvalue, , 11], H.ARL.L2.m[ind.L2.tilt12.hasvalue, , 12], H.ARL.L2.m[ind.L2.tilt13.hasvalue, , 13], H.ARL.L2.m[ind.L2.tilt14.hasvalue, , 14])
  ##  H.ARL.subtop is height above radar level of the top of the storm, equivalent to ECHO TOP HEIGHT,  SINGLE VALUE FOR RADAR DOMAIN
  #H.ARL.subtop                    <- max(H.ARL.max.fullvolume.L2.array.m)     # [m]
  
  # MEHS ('max-estimated-hail-size', Witt et al., 1998) calculation
  #   uses weighting functions in REFL with TEMP-based height thresholds to account for transition zones between rain/hail
  #   defined as layer between 0 and -20 deg [C] isotherm to only include optimal hail growth zone where largest hail likely
  
  ## weighting function W.subT.ofH, from nearby sounding or NWP output
  #W.subT.ofH                   <- array(data = 0L, dim = c(dim(nexrad.L2.tilt1.Z.new)[1], dim(nexrad.L2.tilt1.Z.new)[2], length(theta.radar.volume.L3.deg)))
  #for (j in 1:length(theta.radar.volume.L3.deg)) {
  ##for (j in 4:4) {
  #  #print(j)
  #  for (i in 1:length(nexrad.L3.range.array.m)) {
  #    #print(i)
  #    if (H.ARL.L3.m[i,j] <= H.ARL.sub0) {
  #      W.subT.ofH[i,,j] <- 0
  #    } else if (H.ARL.L3.m[i,j] >= H.ARL.subm20) {
  #      W.subT.ofH[i,,j] <- 1
  #    } else if (H.ARL.L3.m[i,j] > H.ARL.sub0 & H.ARL.L3.m[i,j] < H.ARL.subm20) {
  #      W.subT.ofH[i,,j] <- (H.ARL.L3.m[i,j] - H.ARL.sub0) / (H.ARL.subm20 - H.ARL.sub0)
  #    }  # end of if (H.ARL.m[i,j] <= H.ARL.sub0)
  #    
  #  }    # end of for (i in 1:length(nexrad.L3.range.array.m))
  #  
  #}      # end of for (j in 1:length(theta.radar.volume.deg))
  
  ## weighting function W.subZ, define a transition zone between rain and hail reflectivities
  #W.subZ                     <- array(data = 0L, dim=c(dim(nexrad.L2.tilt1.Z.new)[1], dim(nexrad.L2.tilt1.Z.new)[2], length(theta.radar.volume.L3.deg)))
  #for (j in 1:dim(nexrad.L2.tilt1.Z.new)[2]) {
  #  for (i in 1:dim(nexrad.L2.tilt1.Z.new)[1]) {
  #    if (nexrad.L2.tilt1.Z.new[i,j] <= Z.subL | is.na(nexrad.L2.tilt1.Z.new[i,j])) {
  #      W.subZ[i,j,1] <- 0
  #    } else if (nexrad.L2.tilt1.Z.new[i,j] >= Z.subU) {
  #      W.subZ[i,j,1] <- 1
  #    } else if (nexrad.L2.tilt1.Z.new[i,j] > Z.subL & nexrad.L2.tilt1.Z.new[i,j] < Z.subU) {
  #      W.subZ[i,j,1] <- (nexrad.L2.tilt1.Z.new[i,j] - Z.subL)/(Z.subU - Z.subL)
  #    }  # end of if (nexrad.L2.tilt1.Z.new[i,2] <= Z.subL)
  #  }
  #}
  
  #for (j in 1:dim(nexrad.L2.tilt1.Z.new)[2]) {
  #  for (i in 1:dim(nexrad.L2.tilt2.Z.new)[1]) {   
  #    if (nexrad.L2.tilt2.Z.new[i,j] <= Z.subL | is.na(nexrad.L2.tilt2.Z.new[i,j])) {
  #      W.subZ[i,j,2] <- 0
  #    } else if (nexrad.L2.tilt2.Z.new[i,j] >= Z.subU) {
  #      W.subZ[i,j,2] <- 1
  #    } else if (nexrad.L2.tilt2.Z.new[i,j] > Z.subL & nexrad.L2.tilt2.Z.new[i,j] < Z.subU) {
  #      W.subZ[i,j,2] <- (nexrad.L2.tilt2.Z.new[i,j] - Z.subL)/(Z.subU - Z.subL)
  #    }  # end of if (nexrad.L2.tilt2.Z.new[i,2] <= Z.subL)
  #  }
  #}
  
  #for (j in 1:dim(nexrad.L2.tilt1.Z.new)[2]) {
  #  for (i in 1:dim(nexrad.L2.tilt3.Z.new)[1]) {
  #    if (nexrad.L2.tilt3.Z.new[i,j] <= Z.subL | is.na(nexrad.L2.tilt3.Z.new[i,j])) {
  #      W.subZ[i,j,3] <- 0
  #    } else if (nexrad.L2.tilt3.Z.new[i,j] >= Z.subU) {
  #      W.subZ[i,j,3] <- 1
  #    } else if (nexrad.L2.tilt3.Z.new[i,j] > Z.subL & nexrad.L2.tilt3.Z.new[i,j] < Z.subU) {
  #      W.subZ[i,j,3] <- (nexrad.L2.tilt3.Z.new[i,j] - Z.subL)/(Z.subU - Z.subL)
  #    }  # end of if (nexrad.L2.tilt3.Z.new[i,2] <= Z.subL)
  #  }
  #}
  
  #for (j in 1:dim(nexrad.L2.tilt4.Z.new)[2]) {
  #  for (i in 1:dim(nexrad.L2.tilt4.Z.new)[1]) {  
  #    if (nexrad.L2.tilt4.Z.new[i,j] <= Z.subL | is.na(nexrad.L2.tilt4.Z.new[i,j])) {
  #      W.subZ[i,j*2,4] <- 0
  #    } else if (nexrad.L2.tilt4.Z.new[i,j] >= Z.subU) {
  #      W.subZ[i,j*2,4] <- 1
  #    } else if (nexrad.L2.tilt4.Z.new[i,j] > Z.subL & nexrad.L2.tilt4.Z.new[i,j] < Z.subU) {
  #      W.subZ[i,j*2,4] <- (nexrad.L2.tilt4.Z.new[i,j] - Z.subL)/(Z.subU - Z.subL)
  #    }  # end of if (nexrad.L2.tilt4.Z.new[i,2] <= Z.subL)
  #    W.subZ[i,(j*2)-1,4] <- W.subZ[i,j*2,4]
  #  }
  #}
  
  #for (j in 1:dim(nexrad.L2.tilt5.Z.new)[2]) {
  #  for (i in 1: dim(nexrad.L2.tilt5.Z.new)[1]) {
  #    if (nexrad.L2.tilt5.Z.new[i,j] <= Z.subL | is.na(nexrad.L2.tilt5.Z.new[i,j])) {
  #      W.subZ[i,j*2,5] <- 0
  #    } else if (nexrad.L2.tilt5.Z.new[i,j] >= Z.subU) {
  #      W.subZ[i,j*2,5] <- 1
  #    } else if (nexrad.L2.tilt5.Z.new[i,j] > Z.subL & nexrad.L2.tilt5.Z.new[i,j] < Z.subU) {
  #      W.subZ[i,j*2,5] <- (nexrad.L2.tilt5.Z.new[i,j] - Z.subL)/(Z.subU - Z.subL)
  #    }  # end of if (nexrad.L2.tilt5.Z.new[i,2] <= Z.subL)
  #    W.subZ[i,(j*2)-1,5] <- W.subZ[i,j*2,5]
  #  }
  #}
  
  #for (j in 1:dim(nexrad.L2.tilt6.Z.new)[2]) {
  #  for (i in 1:dim(nexrad.L2.tilt6.Z.new)[1]) {    
  #    if (nexrad.L2.tilt6.Z.new[i,j] <= Z.subL | is.na(nexrad.L2.tilt6.Z.new[i,j])) {
  #      W.subZ[i,j,6] <- 0
  #    } else if (nexrad.L2.tilt6.Z.new[i,j] >= Z.subU) {
  #      W.subZ[i,j,6] <- 1
  #    } else if (nexrad.L2.tilt6.Z.new[i,j] > Z.subL & nexrad.L2.tilt6.Z.new[i,j] < Z.subU) {
  #      W.subZ[i,j,6] <- (nexrad.L2.tilt6.Z.new[i,j] - Z.subL)/(Z.subU - Z.subL)
  #    }  # end of if (nexrad.L2.tilt6.Z.new[i,2] <= Z.subL)
  #    W.subZ[i,(j*2)-1,6] <- W.subZ[i,j*2,6]
  #  }    # end of for (i in 1:length(nexrad.L3.range.array.m))
  #}      # end of for (j in 1:length(theta.radar.volume.array.m))
  
  # THIS COMMENTED STUFF WAS AN ATTEMPT TO UNDERSTAND BELOW ENERGY AND SEV>HAIL>IND COMPUTATIONS
  #ind.H.ARL.betw.sub0.subtop  <- which(H.ARL.L2.m[1:1200,,1:6] >= H.ARL.sub0 & H.ARL.L2.m[1:1200,,1:6] <= H.ARL.subtop) 
