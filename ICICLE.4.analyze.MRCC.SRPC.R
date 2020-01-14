#--------------------------------------------
#
# Name:    ICICLE.4.analyze.MRCC.SRPC.R
#
# Purpose: 1) load mosaic radar (MR) and single radar (SR) RadIA
#          2) make plots and stats
#
# How to run:
#          1) run ICICLE.1.related.data.frames.R to get icicle.df
#          2)
#          3)
#
# Created: 8.20.2019 dserke
#
# Things to do:
# 1. change colorscale on 2x2 INT plots.  Maybe image.plot?
# 2. convert more flights from mdv to nc and cp to pens archive (High ranked flight nums 26, 21, 16)
# 3. test output of loop over each flight.track.df point for accuracy
# 4. is output from F24 SRPC accurate? should return max int in area? explore profiles of int?
# 5. change top row of 2x2 plot to reflect radar images from the mean time of the case period
# 6. setup/plot/output hh:mm on closest algo int value for each flight
# 7. setup scatterplot of FZDZ MCCC vs. SRPC and plot for flight and sumflights
# 8. plot profile of SRPC vs MCCC algos at closest missed approach times/locations
# 9.
#----------------------------------------------

# add libraries
#  data processing
library(foreign)
library(dplyr)
library(magrittr)
library(tidyr)
library(gridExtra)
library(fields)
library(stringr)

# data plotting
library(ggplot2)

#  spatial
library(raster)
library(rasterVis)
library(rgdal)
library(dismo)

#  mapping
library(lattice)
library(latticeExtra)
library(data.table)
#library(tmap)
library(maps)
library(mapdata)
library(maptools)
library(ggmap)

# data input format
library(ncdf4)

# unit conversion
library(NISTunits)

# date manipulation
library(lubridate)

#--------------------
# MAP SERVICE ACQUISITION
#--------------------
#   send Dave's google API account info
register_google(key = "AIzaSyBoSkGvKYWufZCQ33YcxzVfUKA89F55iw4")               # copied directly from Google Console via 'copy' button
#--------------------

#----------------------------------------------
# constants and parameters
#----------------------------------------------
options(max.print=10000) 

# SRPC file input
ind.SRPC.alt                       <- 1                                       # this number of tilts, third index (ex: [1832, 720, 9])
range.to.first.range.gate.km       <- 1.000                                   # [km]
radar.azimuth.resolution.0.5.deg   <- 0.5                                     # every 0.5 deg azim at the lowest tilts

# NSSL mosaic MCCC lat/lon grid bounds, predefined
ICICLE.max.lon                     <- 95                                      # deg W
ICICLE.min.lon                     <- 83                                      # deg W
ICICLE.max.lat                     <- 48                                      # deg N
ICICLE.min.lat                     <- 36                                      # deg N

# MCCC file input 
ind.nssl.mosaic.alt                <- 3                                       # 3rd alt in array below
MCCC.alt.array                     <- c(0.5, 0.75, 1.00, 1.25, 1.50, 1.75, 2.00, 2.25, 2.50, 2.75, 3.00, 3.50, 4.00, 4.50, 5.00, 5.50, 6.00, 6.50, 7.00, 7.50, 8.00, 8.50, 9.00, 10.00)
MCCC.lat.array                     <- seq(from =  ICICLE.min.lat, to =  ICICLE.max.lat-0.01, by = 0.01)
MCCC.lon.array                     <- seq(from = -ICICLE.max.lon, to = -ICICLE.min.lon-0.01, by = 0.01)

# indices to plot RadIA MCCC and SRPC profiles
#   these were chosen for the 20190227 2:10Z case as they fall in the middle of the return pattern
ind.lon                            <- 360
ind.lat                            <- 600

radius.earth.km                    <- 6378.1                                  # [km]
degPERcircle                       <- 360                                     # [deg]
kmPERm                             <- 0.001                                   # convert m to km
kmPERrangegate                     <- 0.250                                   #

#----------------------------------------------
# define data paths,list data file in path
#----------------------------------------------
# define NEXRAD site location csv file path
nexrad.site.dataframetxt.dir       <- file.path("/d1/serke/projects/case_studies/SNOWIE/data/RadIA_data/nexrad_site_data/")
# define NRC CNV flight track (planet) csv file path
flight.track.planet.dir            <- file.path("/d1/serke/projects/case_studies/ICICLE_2018/data/flight_track_planet/")
flight.track.planet.listfiles      <- list.files(path = flight.track.planet.dir, include.dirs = FALSE)
# define output array file path
closest.int.val.array.dir          <- file.path("/d1/serke/projects/case_studies/ICICLE_2018/data/R_processing_output/")

#----------------------------------------------
# load the NEXRAD site location text file
#----------------------------------------------
# load csv file
NEXRAD.site.df                     <- read.csv(paste(nexrad.site.dataframetxt.dir, "nexrad_site.csv", sep = ""), header = FALSE, sep = ",", dec = ".", stringsAsFactors=FALSE)
colnames(NEXRAD.site.df)           <- c("NCDCID", "ICAO", "WBAN", "radname", "COUNTRY", "STATE", "COUNTY", "lat", "lon", "elev", "GMTdiff", "STN_TYPE")
head(NEXRAD.site.df)

#----------------------------------------------
# loop over length of icicle.df
#----------------------------------------------
ii                                 <- 0                      # a counter 
fl.num.wdata.avail                 <- c(7, 17, 21, 24)       # these fl nums have nc data available
ind.fl.num.wdata.avail             <- array(data = NA, dim = length(fl.num.wdata.avail))
for (p in 1:length(fl.num.wdata.avail)) {
  ind.fl.num.wdata.avail[p] <- which(icicle.df$fl.num == fl.num.wdata.avail[p])
}
print(ind.fl.num.wdata.avail)                                # the index nums in icicle.df of all fl nums with available nc data

# loop over only the indices of fl nums with available nc data files
#for (i in ind.fl.num.wdata.avail) { # for all available flights
for (i in 10:10) {
    
  ii                                 <- ii + 1
  yyyymmdd.for.case                  <- paste("", icicle.df$date.posix[i], sep="")
  
  #----------------------------------------------
  # define NEXRAD radar-specific information for case i
  #----------------------------------------------
  #   for primary (pri) radar
  radar.for.case                       <- paste("", icicle.df$pri.radar.name[i], sep="")
  ind.radar                            <- which(NEXRAD.site.df$ICAO == paste(" ", icicle.df$pri.radar.name[i], sep=""))
  radar.lat.deg                        <- NEXRAD.site.df$lat[ind.radar]
  radar.lon.deg                        <- NEXRAD.site.df$lon[ind.radar]
  radar.alt.MSL.ft                     <- NEXRAD.site.df$elev[ind.radar]
  #   for secondary (sec) radar
  #   NOTE: NO PLOTTING/PROCESSING IS YET SET UP FOR SECONDARY RADAR FOR A GIVEN CASE.  NEED TO CONSIDER THIS
  if (!is.na(paste("", icicle.df$sec.radar.name[i], sep=""))) {
    radar.sec.for.case                 <- paste("", icicle.df$sec.radar.name[i], sep="")
    ind.radar.sec                      <- which(NEXRAD.site.df$ICAO == paste(" ", icicle.df$sec.radar.name[i], sep=""))
    radar.sec.lat.deg                  <- NEXRAD.site.df$lat[ind.radar.sec]
    radar.sec.lon.deg                  <- NEXRAD.site.df$lon[ind.radar.sec]
    radar.sec.alt.MSL.ft               <- NEXRAD.site.df$elev[ind.radar.sec]
  } # end of if (!is.na()...)
  
  #----------------------------------------------
  # define radar data paths for case yyyymmdd.for.case then list data file in path
  #----------------------------------------------
  # define RadIA SRPC nc interest directories
  radia.SRPC.nc.ALL.dir              <- file.path(paste("/d1/serke/projects/case_studies/ICICLE_2018/data/nc_radia_SRPC/", radar.for.case, "/polar/", yyyymmdd.for.case, "/", sep=""))
  radia.SRPC.nc.ALL.listfiles        <- list.files(path = radia.SRPC.nc.ALL.dir , include.dirs = FALSE) 
  # define RadIA MRCC nc interest directories
  radia.MCCC.nc.FZDZ.dir             <- file.path(paste("/d1/serke/projects/case_studies/ICICLE_2018/data/nc_radia_MCCC/FZDZ/",   yyyymmdd.for.case, "/", sep=""))
  radia.MCCC.nc.FZDZ.listfiles       <- list.files(path = radia.MCCC.nc.FZDZ.dir , include.dirs = FALSE) 
  radia.MCCC.nc.SLW.dir              <- file.path(paste("/d1/serke/projects/case_studies/ICICLE_2018/data/nc_radia_MCCC/SLW/",    yyyymmdd.for.case, "/", sep=""))
  radia.MCCC.nc.SLW.listfiles        <- list.files(path = radia.MCCC.nc.SLW.dir , include.dirs = FALSE)
  radia.MCCC.nc.MIXPHA.dir           <- file.path(paste("/d1/serke/projects/case_studies/ICICLE_2018/data/nc_radia_MCCC/MIXPHA/", yyyymmdd.for.case, "/", sep=""))
  radia.MCCC.nc.MIXPHA.listfiles     <- list.files(path = radia.MCCC.nc.MIXPHA.dir , include.dirs = FALSE)
  
  #----------------------------------------------
  # load the flight track csv file from planet and calc mean lat and lon for later plotting
  #----------------------------------------------
  ind.flttrk.file                    <- 0
  for (k in 1:length(flight.track.planet.listfiles)) {
    torf.flttrk.file <- str_detect(flight.track.planet.listfiles[k], paste("F", icicle.df$fl.num[i], "_NRC", sep=""))
    if (torf.flttrk.file == TRUE) {ind.flttrk.file <- k}
  }  # end of for ( k in 1:length())
  flight.track.df                    <- read.csv(paste(flight.track.planet.dir, flight.track.planet.listfiles[ind.flttrk.file], sep = ""), header = FALSE, sep = ",", dec = ".", stringsAsFactors=FALSE)
  colnames(flight.track.df)          <- c("type", "ac.date.posix", "lat", "lon", "alt.m", "", "alt.ft", "", "radar.alt.ft", "ground.spd.mps", "true.spd.mps", "ind.spd.kts", "mach.num", "vert.vel.mps", "true.hdg.deg", "track.deg", "drift.deg", "pitch.deg", "roll.deg", "side.slip.deg", "angle.atk.deg", "amb.temp.deg.C", "dewpt.temp.deg.C", "", "", "", "?5", "?6")
  head(flight.track.df)  
  lat.flttrk.mean                    <- mean(flight.track.df$lat, na.rm = TRUE)
  lon.flttrk.mean                    <- mean(flight.track.df$lon, na.rm = TRUE)
  
  # use pri.radar.srt and pri.radar.end hh:mm to loop through available MCCC and SRPC nc files
  posix.hhmm.srt                     <- as.POSIXct(paste(icicle.df$date.posix[i], " ", icicle.df$pri.radar.srt[i], sep=""))
  posix.hhmm.end                     <- as.POSIXct(paste(icicle.df$date.posix[i], " ", icicle.df$pri.radar.end[i], sep=""))
  if (posix.hhmm.srt > posix.hhmm.end) {
    posix.hhmm.end  <- posix.hhmm.end + (24*60*60)  # add 1 day to the posix time
  }  # end of if (posix.hhmm.srt > ...)
  diff.time.hr                       <- round(difftime(posix.hhmm.end, posix.hhmm.srt))
  
  #---------------------------------------------
  # calculate (hh*60)*mm for radar nc files
  #   use this value for MCCC, SRPC to find closest time to each flight track time
  #---------------------------------------------
  #   calc tot mins of each SRPC file times (hh*mm)
  mm.tot.SRPC.file                   <- array(data = NA, dim = c(length(radia.SRPC.nc.ALL.listfiles), 1))
  for (a in 1:length(radia.SRPC.nc.ALL.listfiles)) { 
    hh.SRPC.file             <- as.numeric(substr(radia.SRPC.nc.ALL.listfiles[a], 10, 11))
    mm.SRPC.file             <- as.numeric(substr(radia.SRPC.nc.ALL.listfiles[a], 12, 13))
    mm.tot.SRPC.file[a]      <- (hh.SRPC.file * 60) + mm.SRPC.file
  }  # end of for (a...)
  
  mm.tot.MCCC.FZDZ.file              <- array(data = NA, dim = c(length(radia.MCCC.nc.FZDZ.listfiles), 1))
  #   calc tot mins of each MCCC FZDZ file times (hh*mm)
  for (b in 1:length(radia.MCCC.nc.FZDZ.listfiles)) {
    hh.MCCC.file             <- as.numeric(substr(radia.MCCC.nc.FZDZ.listfiles[b], 10, 11))
    mm.MCCC.file             <- as.numeric(substr(radia.MCCC.nc.FZDZ.listfiles[b], 12, 13))
    mm.tot.MCCC.FZDZ.file[b] <- (hh.MCCC.file * 60) + mm.MCCC.file
  }  # end of for (b...)
  
  mm.tot.MCCC.SLW.file               <- array(data = NA, dim = c(length(radia.MCCC.nc.SLW.listfiles), 1))
  #   calc tot mins of each MCCC FZDZ file times (hh*mm)
  for (c in 1:length(radia.MCCC.nc.FZDZ.listfiles)) {
    hh.MCCC.file             <- as.numeric(substr(radia.MCCC.nc.SLW.listfiles[c], 10, 11))
    mm.MCCC.file             <- as.numeric(substr(radia.MCCC.nc.SLW.listfiles[c], 12, 13))
    mm.tot.MCCC.SLW.file[c]  <- (hh.MCCC.file * 60) + mm.MCCC.file
  }  # end of for (c...)
  
  mm.tot.MCCC.MIXPHA.file            <- array(data = NA, dim = c(length(radia.MCCC.nc.MIXPHA.listfiles), 1))
  #   calc tot mins of each MCCC FZDZ file times (hh*mm)
  for (d in 1:length(radia.MCCC.nc.MIXPHA.listfiles)) {
    hh.MCCC.file             <- as.numeric(substr(radia.MCCC.nc.MIXPHA.listfiles[d], 10, 11))
    mm.MCCC.file             <- as.numeric(substr(radia.MCCC.nc.MIXPHA.listfiles[d], 12, 13))
    mm.tot.MCCC.MIXPHA.file[d]<- (hh.MCCC.file * 60) + mm.MCCC.file
  }  # end of for (d...)
  
  #browser()
  
  #----------------------------------------------
  # print some info on the loaded flight fields to the GUI
  #----------------------------------------------
  print(paste("Flight #", fl.num.wdata.avail[ii], ", date = ", yyyymmdd.for.case, sep=""))
  print(paste("  Primary radar = ",                            radar.for.case, sep=""))
  print(paste("    Srt time in domain = ",                     posix.hhmm.srt, sep=""))
  print(paste("    End time in domain = ",                     posix.hhmm.end, sep=""))
  print(paste("    Num hrs in radar domain = ",                diff.time.hr))
  print(paste("    Num SRPC ALL files avail for flight = ",    length(radia.SRPC.nc.ALL.listfiles)))
  print(paste("    Num MCCC FZDZ files avail for flight = ",   length(radia.MCCC.nc.FZDZ.listfiles)))
  print(paste("    Num MCCC SLW files avail for flight = ",    length(radia.MCCC.nc.SLW.listfiles)))
  print(paste("    Num MCCC MIXPHA files avail for flight = ", length(radia.MCCC.nc.MIXPHA.listfiles)))
  print(paste("  Flight track file = ",                        flight.track.planet.listfiles[ind.flttrk.file], sep=""))
  print(paste("    Num pts in file = ",                        dim(flight.track.df)[1]))
  
  #-----------------------------------------------
  # loop through all flight track times and find closest SRPC and MCCC data file to given flight time
  #-----------------------------------------------
  ind.SRPC.file.closest.last        <- 0  # start with no last index
  ind.MCCC.FZDZ.file.closest.last   <- 0
  ind.MCCC.SLW.file.closest.last    <- 0
  ind.MCCC.MIXPHA.file.closest.last <- 0
  
  #################################################
  # MUST EDIT LINES BELOW FOR EACH DIFFERENT FLIGHT
  #################################################
  #ind.flttrk.times                  <- 409:588  # for F7
  #ind.flttrk.times                  <- 31:331   # for F17
  #ind.flttrk.times                  <- 256:436  # for F21
  #ind.flttrk.times                  <- 620:668  # for F24, 0207-0238 Z
  ind.flttrk.times                  <- 493
  
  # initialize arrays for each field, one entry per value in flight track 
  length.flight.track               <- dim(flight.track.df)[1]
  ind.closest.MCCC.lat              <- array(data = 0L, dim = c(dim(flight.track.df)[1]))
  ind.closest.MCCC.lon              <- array(data = 0L, dim = c(dim(flight.track.df)[1]))
  ind.closest.MCCC.alt              <- array(data = 0L, dim = c(dim(flight.track.df)[1]))
  FZDZ.closest.MCCC                 <- array(data = NA, dim = c(dim(flight.track.df)[1]))
  SLW.closest.MCCC                  <- array(data = NA, dim = c(dim(flight.track.df)[1]))
  MIXPHA.closest.MCCC               <- array(data = NA, dim = c(dim(flight.track.df)[1]))
  X.rad                             <- array(data = 0L, dim = c(dim(flight.track.df)[1]))
  Y.rad                             <- array(data = 0L, dim = c(dim(flight.track.df)[1]))
  bearing.radartoac.deg             <- array(data = 0L, dim = c(dim(flight.track.df)[1]))
  bearing.radartoac.ind             <- array(data = 0L, dim = c(dim(flight.track.df)[1]))
  a                                 <- array(data = 0L, dim = c(dim(flight.track.df)[1]))
  c                                 <- array(data = 0L, dim = c(dim(flight.track.df)[1]))
  dist.radartoac.km                 <- array(data = 0L, dim = c(dim(flight.track.df)[1]))
  alpha.deg                         <- array(data = 0L, dim = c(dim(flight.track.df)[1]))
  range.radartoac.km                <- array(data = 0L, dim = c(dim(flight.track.df)[1]))
  range.radartoac.ind               <- array(data = 0L, dim = c(dim(flight.track.df)[1]))
  theta.radartoac.deg               <- array(data = 0L, dim = c(dim(flight.track.df)[1]))
  theta.radartoac.ind               <- array(data = 0L, dim = c(dim(flight.track.df)[1]))
  FZDZ.closest.SRPC                 <- array(data = NA, dim = c(dim(flight.track.df)[1]))
  SLW.closest.SRPC                  <- array(data = NA, dim = c(dim(flight.track.df)[1]))
  MIXPHA.closest.SRPC               <- array(data = NA, dim = c(dim(flight.track.df)[1]))
  hhmm.closest                      <- array(data = NA, dim = c(length(ind.flttrk.times)))
    
  for (j in ind.flttrk.times) {
    
    # calc tot mins of flight track time (hh*mm)
    #print(j)
    hh.flight.track              <- hour(flight.track.df$ac.date.posix[j])
    mm.flight.track              <- minute(flight.track.df$ac.date.posix[j])
    hhmm.closest[j]              <- paste(hh.flight.track, ":", mm.flight.track, sep="")
    mm.tot.flight.track          <- (hh.flight.track * 60) + mm.flight.track
    
    # find the closest hh*mm time for SRPC and MCCC
    ind.SRPC.file.closest        <- which.min(abs(mm.tot.flight.track - mm.tot.SRPC.file))
    ind.MCCC.FZDZ.file.closest   <- which.min(abs(mm.tot.flight.track - mm.tot.MCCC.FZDZ.file))
    ind.MCCC.SLW.file.closest    <- which.min(abs(mm.tot.flight.track - mm.tot.MCCC.SLW.file))
    ind.MCCC.MIXPHA.file.closest <- which.min(abs(mm.tot.flight.track - mm.tot.MCCC.MIXPHA.file))
    
    #----------------------------------------------
    # load RadIA interest files, if not already loaded
    #----------------------------------------------
    # SRPC file containing FZDZ, SLW, MIXPHA and PLATES
    if (ind.SRPC.file.closest != ind.SRPC.file.closest.last) {
      radia.SRPC.nc.ALL.filename        <- paste(radia.SRPC.nc.ALL.dir, radia.SRPC.nc.ALL.listfiles[ind.SRPC.file.closest], sep="")
      print(paste("loading ", radia.SRPC.nc.ALL.filename, sep=""))
      if (substr(radia.SRPC.nc.ALL.filename, nchar(radia.SRPC.nc.ALL.filename)-1, nchar(radia.SRPC.nc.ALL.filename)) == "gz") {
        print("File is gzipped. Unzipping...")
        untar(radia.SRPC.nc.ALL.filename)
        radia.SRPC.nc.ALL.filename <- substr(radia.SRPC.nc.ALL.filename, 1, nchar(radia.SRPC.nc.ALL.filename)-3)
      }  # end of if (substr(radia.SRPC.nc.ALL.filename, ...))
      radia.SRPC.nc.ALL                 <- nc_open(radia.SRPC.nc.ALL.filename, write = FALSE, verbose = FALSE)
      theta.deg                         <- radia.SRPC.nc.ALL$dim[2]$height$vals
      #print(paste("The file has", radia.SRPC.nc.ALL$nvars, "variables"))
      radia.SRPC.ALL.var.num           <- seq(1, radia.SRPC.nc.ALL$nvars, by=1)
      for (s in 1:length(radia.SRPC.ALL.var.num)) {
        radia.SRPC.ALL.nam <- paste("v", radia.SRPC.ALL.var.num[s], sep = "")
        assign(radia.SRPC.ALL.nam, radia.SRPC.nc.ALL$var[[radia.SRPC.ALL.var.num[s]]])
      }  # end of for (i in 1:length(radia.SRPC.ALL.var.num))
      if (radia.SRPC.nc.ALL$nvars == 11) {
        FZDZ.SRPC                                <- ncvar_get( radia.SRPC.nc.ALL, v7 )
        SLW.SRPC                                 <- ncvar_get( radia.SRPC.nc.ALL, v8 )
        MIXPHA.SRPC                              <- ncvar_get( radia.SRPC.nc.ALL, v9 )
        PLATES.SRPC                              <- ncvar_get( radia.SRPC.nc.ALL, v10 )
      } else if (radia.SRPC.nc.ALL$nvars == 1) {
        FZDZ.SRPC                                <- ncvar_get( radia.SRPC.nc.ALL, v1 )
      }  # end of if (radia.SRPC.nc.ALL$nvars == 5)
      #print(paste("V1 has name", v1$name))
      nc_close(radia.SRPC.nc.ALL)
    }    # end of if (ind.SRPC.file.closest ...)

    # MCCC file containing FZDZ
    if (ind.MCCC.FZDZ.file.closest != ind.MCCC.FZDZ.file.closest.last) {  
      radia.MCCC.nc.FZDZ.filename       <- paste(radia.MCCC.nc.FZDZ.dir, radia.MCCC.nc.FZDZ.listfiles[ind.MCCC.FZDZ.file.closest], sep="")
      print(paste("loading ", radia.MCCC.nc.FZDZ.filename, sep=""))
      if (substr(radia.MCCC.nc.FZDZ.filename, nchar(radia.MCCC.nc.FZDZ.filename)-1, nchar(radia.MCCC.nc.FZDZ.filename)) == "gz") {
        print("File is gzipped. Unzipping...")
        untar(radia.MCCC.nc.FZDZ.filename)
        radia.MCCC.nc.FZDZ.filename <- substr(radia.MCCC.nc.FZDZ.filename, 1, nchar(radia.MCCC.nc.FZDZ.filename)-3)
      }  # end of if ()...
      radia.MCCC.nc.FZDZ                <- nc_open(radia.MCCC.nc.FZDZ.filename, write = FALSE, verbose = FALSE)
      #print(paste("The file has", radia.MCCC.nc.FZDZ$nvars, "variables"))
      radia.MCCC.FZDZ.var.num                <- seq(1, radia.MCCC.nc.FZDZ$nvars, by=1)
      for (s in 1:length(radia.MCCC.FZDZ.var.num)) {
        radia.MCCC.FZDZ.nam <- paste("v", radia.MCCC.FZDZ.var.num[s], sep = "")
        assign(radia.MCCC.FZDZ.nam, radia.MCCC.nc.FZDZ$var[[radia.MCCC.FZDZ.var.num[s]]])
      }  # end of for ()...
      if (radia.MCCC.nc.FZDZ$nvars >= 9) {
        FZDZ.MCCC                                <- ncvar_get( radia.MCCC.nc.FZDZ, v6 )
        meanDBZ.MCCC                             <- ncvar_get( radia.MCCC.nc.FZDZ, v7 )
        sdevDBZ.MCCC                             <- ncvar_get( radia.MCCC.nc.FZDZ, v8 )
        TDBZ.MCCC                                <- ncvar_get( radia.MCCC.nc.FZDZ, v9 )
      } else if (radia.MCCC.nc.FZDZ$nvars == 1) {
        FZDZ.MCCC                                <- ncvar_get( radia.MCCC.nc.FZDZ, v1 )
      }  # end of if ()...
      #print(paste("V1 has name", v1$name))
      nc_close(radia.MCCC.nc.FZDZ)
      FZDZ.MCCC.df <- data.frame(FZDZ.MCCC[ , , ind.nssl.mosaic.alt])
      #head(FZDZ.MCCC.df)
    }    # end of if (ind.MCCC.file.closest...)

    # MCCC file containing SLW
    if (ind.MCCC.SLW.file.closest != ind.MCCC.SLW.file.closest.last) {
      radia.MCCC.nc.SLW.filename       <- paste(radia.MCCC.nc.SLW.dir, radia.MCCC.nc.SLW.listfiles[ind.MCCC.SLW.file.closest], sep="")
      print(paste("loading ", radia.MCCC.nc.SLW.filename, sep=""))
      if (substr(radia.MCCC.nc.SLW.filename, nchar(radia.MCCC.nc.SLW.filename)-1, nchar(radia.MCCC.nc.SLW.filename)) == "gz") {
        print("File is gzipped. Unzipping...")
        untar(radia.MCCC.nc.SLW.filename)
        radia.MCCC.nc.SLW.filename <- substr(radia.MCCC.nc.SLW.filename, 1, nchar(radia.MCCC.nc.SLW.filename)-3)
      }  # end of if () ...
      radia.MCCC.nc.SLW                     <- nc_open(radia.MCCC.nc.SLW.filename, write = FALSE, verbose = FALSE)
      #print(paste("The file has", radia.MCCC.nc.SLW$nvars, "variables"))
      radia.MCCC.SLW.var.num                <- seq(1, radia.MCCC.nc.SLW$nvars, by=1)
      for (s in 1:length(radia.MCCC.SLW.var.num)) {
        radia.MCCC.SLW.nam <- paste("v", radia.MCCC.SLW.var.num[s], sep = "")
        assign(radia.MCCC.SLW.nam, radia.MCCC.nc.SLW$var[[radia.MCCC.SLW.var.num[s]]])
      }  # end of for () ...
      if (radia.MCCC.nc.SLW$nvars >= 10) {
        SLW.MCCC                                <- ncvar_get( radia.MCCC.nc.SLW, v6 )
        meanZDR.MCCC                            <- ncvar_get( radia.MCCC.nc.SLW, v7 )
        sdevZDR.MCCC                            <- ncvar_get( radia.MCCC.nc.SLW, v8 )
        meanKDP.MCCC                            <- ncvar_get( radia.MCCC.nc.SLW, v9 )
        sdevKDP.MCCC                            <- ncvar_get( radia.MCCC.nc.SLW, v10)
      } else if (radia.MCCC.nc.SLW$nvars == 1) {
        SLW.MCCC                                <- ncvar_get( radia.MCCC.nc.SLW, v1 )
      }  # end of if ()...
      #print(paste("V1 has name", v1$name))
      nc_close(radia.MCCC.nc.SLW)
      SLW.MCCC.df <- data.frame(SLW.MCCC[ , , ind.nssl.mosaic.alt])
    }    # end of if (ind.MCCC.SLW.file.closest ~= ...)

    # MCCC file containing MIXPHA
    if (ind.MCCC.MIXPHA.file.closest != ind.MCCC.MIXPHA.file.closest.last) {
      radia.MCCC.nc.MIXPHA.filename       <- paste(radia.MCCC.nc.MIXPHA.dir, radia.MCCC.nc.MIXPHA.listfiles[ind.MCCC.MIXPHA.file.closest], sep="")
      print(paste("loading ", radia.MCCC.nc.MIXPHA.filename, sep=""))
      if (substr(radia.MCCC.nc.MIXPHA.filename, nchar(radia.MCCC.nc.MIXPHA.filename)-1, nchar(radia.MCCC.nc.MIXPHA.filename)) == "gz") {
        print("File is gzipped. Unzipping...")
        untar(radia.MCCC.nc.MIXPHA.filename)
        radia.MCCC.nc.MIXPHA.filename <- substr(radia.MCCC.nc.MIXPHA.filename, 1, nchar(radia.MCCC.nc.MIXPHA.filename)-3)
      }  # end of if () ...
      radia.MCCC.nc.MIXPHA                <- nc_open(radia.MCCC.nc.MIXPHA.filename, write = FALSE, verbose = FALSE)
      #print(paste("The file has", radia.MCCC.nc.MIXPHA$nvars, "variables"))
      radia.MCCC.MIXPHA.var.num                <- seq(1, radia.MCCC.nc.MIXPHA$nvars, by=1)
      for (s in 1:length(radia.MCCC.MIXPHA.var.num)) {
        radia.MCCC.MIXPHA.nam <- paste("v", radia.MCCC.MIXPHA.var.num[s], sep = "")
        assign(radia.MCCC.MIXPHA.nam, radia.MCCC.nc.MIXPHA$var[[radia.MCCC.MIXPHA.var.num[s]]])
      }  # end of for ()...
      if (radia.MCCC.nc.MIXPHA$nvars >= 9) {
        MIXPHA.MCCC                               <- ncvar_get( radia.MCCC.nc.MIXPHA, v6 )
        meanZDR.MCCC                              <- ncvar_get( radia.MCCC.nc.MIXPHA, v7 )
        meanBDZ.MCCC                              <- ncvar_get( radia.MCCC.nc.MIXPHA, v8 )
        TEMP_MIXPHA.MCCC                          <- ncvar_get( radia.MCCC.nc.MIXPHA, v9 )
      } else if (radia.MCCC.nc.MIXPHA$nvars == 1) {
        MIXPHA.MCCC                               <- ncvar_get( radia.MCCC.nc.MIXPHA, v1 )
      }  # end of if ()...
      #print(paste("V1 has name", v1$name))
      nc_close(radia.MCCC.nc.MIXPHA)
      MIXPHA.MCCC.df <- data.frame(MIXPHA.MCCC[ , , ind.nssl.mosaic.alt])
    }    # end of if (ind.MCCC.MIXPHA.file.closest ~= ...)
  
    # define indices of the last closest file, to see if a new nc file needs to be loaded or not
    ind.SRPC.file.closest.last        <- ind.SRPC.file.closest
    ind.MCCC.FZDZ.file.closest.last   <- ind.MCCC.FZDZ.file.closest
    ind.MCCC.SLW.file.closest.last    <- ind.MCCC.SLW.file.closest
    ind.MCCC.MIXPHA.file.closest.last <- ind.MCCC.MIXPHA.file.closest
  
    #-------------------------------------
    # Manipulate the data
    #-------------------------------------
    # find indices of closest MCCC lat/lon to each AC lat/lon
    ind.closest.MCCC.lat[j]           <- which.min(abs(flight.track.df$lat[j]   - MCCC.lat.array))
    ind.closest.MCCC.lon[j]           <- which.min(abs(flight.track.df$lon[j]   - MCCC.lon.array))
    ind.closest.MCCC.alt[j]           <- which.min(abs(flight.track.df$alt.m[j] - MCCC.alt.array*1000))
    FZDZ.closest.MCCC[j]              <- FZDZ.MCCC[   ind.closest.MCCC.lon[j], ind.closest.MCCC.lat[j], ind.closest.MCCC.alt[j]]
    SLW.closest.MCCC[j]               <- SLW.MCCC[    ind.closest.MCCC.lon[j], ind.closest.MCCC.lat[j], ind.closest.MCCC.alt[j]]
    MIXPHA.closest.MCCC[j]            <- MIXPHA.MCCC[ ind.closest.MCCC.lon[j], ind.closest.MCCC.lat[j], ind.closest.MCCC.alt[j]]
  
    # find range/azimuth indices of closest SRPC lat/lon to each AC lat/lon
    # CALCULATE BEARING, HORIZ DISTANCE, ALPHA, AND RANGE FROM RADAR TO AC

    # calculate bearing (degrees) between two lat/lon points
    #β = atan2(X,Y) is bearing from X to Y
    #X = cos θb * sin ∆L
    #Y = cos θa * sin θb – sin θa * cos θb * cos ∆L
    X.rad[j]                          <- cos(NISTdegTOradian(flight.track.df$lat[j])) * sin(abs(NISTdegTOradian(radar.lon.deg)-NISTdegTOradian(flight.track.df$lon[j])))
    Y.rad[j]                          <- cos(NISTdegTOradian(radar.lat.deg)) * sin(NISTdegTOradian(flight.track.df$lat[j])) - sin(NISTdegTOradian(radar.lat.deg)) * cos(NISTdegTOradian(flight.track.df$lat[j])) * cos(abs(NISTdegTOradian(radar.lon.deg)-NISTdegTOradian(flight.track.df$lon[j])))
    bearing.radartoac.deg[j]          <- 360 - NISTradianTOdeg(atan2(X.rad[j], Y.rad[j]))
  
    # calculate bearing index of ac to radar
    bearing.radartoac.ind[j] <- round(bearing.radartoac.deg[j] / radar.azimuth.resolution.0.5.deg)
    if (bearing.radartoac.ind[j] == 0) {
      bearing.radartoac.ind[j] <- 720
    }  # end of if (bearing.radartoac.ind[m] )

    #calculate great circle distance along earth's surface between two lat/lon pts
    #a = sin²(Δφ/2) + cos φ1 ⋅ cos φ2 ⋅ sin²(Δλ/2)
    #c = 2 ⋅ atan2( √a, √(1−a) )
    #d = R ⋅ c
    a[j]                        <- sin(abs(NISTdegTOradian(radar.lat.deg) - NISTdegTOradian(flight.track.df$lat[j])) / 2) * sin(abs(NISTdegTOradian(radar.lat.deg)-NISTdegTOradian(flight.track.df$lat[j])) / 2) + (cos(NISTdegTOradian(radar.lat.deg)) * cos(NISTdegTOradian(flight.track.df$lat[j])) * sin(abs(NISTdegTOradian(radar.lon.deg) - NISTdegTOradian(flight.track.df$lon[j])) / 2) * sin(abs(NISTdegTOradian(radar.lon.deg)-NISTdegTOradian(flight.track.df$lon[j])) / 2) )
    c[j]                        <- 2 * atan2(sqrt(a[j]), sqrt(1-a[j]))
    dist.radartoac.km[j]        <- radius.earth.km * c[j]

    # calc alpha.deg, which is angle between radar and ac measured from center of earth
    alpha.deg[j]                <- (dist.radartoac.km[j] * degPERcircle) / (2 * pi * radius.earth.km)
  
    # calculate range.radartoac.km
    range.radartoac.km[j]       <- sqrt((radius.earth.km * radius.earth.km) + (radius.earth.km + (flight.track.df$alt.m[j] * kmPERm))^2 - 2 * radius.earth.km * (radius.earth.km + (flight.track.df$alt.m[j] * kmPERm)) * cos(NISTdegTOradian(alpha.deg[j])))

    # calculate range index of ac from radar
    range.radartoac.ind[j]      <- round((range.radartoac.km[j] - range.to.first.range.gate.km) / kmPERrangegate)
  
    # calculate theta (zenith angle) of ac from radar
    theta.radartoac.deg[j]      <- atan((flight.track.df$alt.m[j] * kmPERm) / dist.radartoac.km[j])
  
    # calculate
    theta.radartoac.ind[j]      <- 1

    # rasterize
    #   convert fields to raster
    FZDZ.MCCC.raster            <- raster(FZDZ.MCCC  [ , , ind.closest.MCCC.alt[j]])   # ind 3 = 1.0 km cart height
    SLW.MCCC.raster             <- raster(SLW.MCCC   [ , , ind.closest.MCCC.alt[j]])
    MIXPHA.MCCC.raster          <- raster(MIXPHA.MCCC[ , , ind.closest.MCCC.alt[j]])
    FZDZ.SRPC.raster            <- raster(FZDZ.SRPC  [ , , 1])                         # ind 1 = first tilt angle in volume = 0.5 degree
    SLW.SRPC.raster             <- raster(SLW.SRPC   [ , , 1])
    MIXPHA.SRPC.raster          <- raster(MIXPHA.SRPC[ , , 1])
    
    #   define x/y extent
    extent(FZDZ.MCCC.raster)    <- extent(ICICLE.min.lat, ICICLE.max.lat, -ICICLE.max.lon, -ICICLE.min.lon)
    extent(SLW.MCCC.raster)     <- extent(ICICLE.min.lat, ICICLE.max.lat, -ICICLE.max.lon, -ICICLE.min.lon)
    extent(MIXPHA.MCCC.raster)  <- extent(ICICLE.min.lat, ICICLE.max.lat, -ICICLE.max.lon, -ICICLE.min.lon)
    extent(FZDZ.SRPC.raster)    <- extent(0.00, 360.00, -458.75, -1)
    extent(SLW.SRPC.raster)     <- extent(0.00, 360.00, -458.75, -1)
    extent(MIXPHA.SRPC.raster)  <- extent(0.00, 360.00, -458.75, -1)
    
    #   set raster CRS values
    crs(FZDZ.MCCC.raster)       <- "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"
    crs(SLW.MCCC.raster)        <- "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"
    crs(MIXPHA.MCCC.raster)     <- "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"
    crs(FZDZ.SRPC.raster)       <- "+init=epsg:3995"
    crs(SLW.SRPC.raster)        <- "+init=epsg:3995"
    crs(MIXPHA.SRPC.raster)     <- "+init=epsg:3995"
    
    #   print some raster values
    #FZDZ.MCCC.raster
    #SLW.MCCC.raster
    #MIXPHA.MCCC.raster
    #FZDZ.SRPC.raster
    #SLW.SRPC.raster
    #MIXPHA.SRPC.raster
    
    # find 
    ind.halfdegs.center.azim    <- -5:5
    ind.rangedate.center.range  <- -10:10
    #if (bearing.radartoac.ind >= dim(FZDZ.SRPC)[2]) {
    #}
    FZDZ.closest.SRPC[j]        <- mean(FZDZ.SRPC  [range.radartoac.ind[j]+ind.halfdegs.center.azim, bearing.radartoac.ind[j]+ind.rangedate.center.range, 1:9], na.rm=TRUE)
    SLW.closest.SRPC[j]         <- mean(SLW.SRPC   [range.radartoac.ind[j]+ind.halfdegs.center.azim, bearing.radartoac.ind[j]+ind.rangedate.center.range, 1:9], na.rm=TRUE)
    MIXPHA.closest.SRPC[j]      <- mean(MIXPHA.SRPC[range.radartoac.ind[j]+ind.halfdegs.center.azim, bearing.radartoac.ind[j]+ind.rangedate.center.range, 1:9], na.rm=TRUE)
    
    ## find closest radar tilt angle to theta.radartoac.deg
    ##ind.radar.closest.tilt <- which.min(abs(volume.tilts-theta.radartoac.deg))
    #FZDZ.closest.SRPC[j]        <- FZDZ.SRPC[   range.radartoac.ind[j], bearing.radartoac.ind[j], theta.radartoac.ind[j]]
    #SLW.closest.SRPC[j]         <- SLW.SRPC[    range.radartoac.ind[j], bearing.radartoac.ind[j], theta.radartoac.ind[j]]
    #MIXPHA.closest.SRPC[j]      <- MIXPHA.SRPC[ range.radartoac.ind[j], bearing.radartoac.ind[j], theta.radartoac.ind[j]]
  
  }  # end of for (j in 1: length(flight.track))

  # set FZDZ.closest.SRPC to NA where AC was not in return
  #ind.FZDZ.closest.MCCC.isNA <- which(is.na(FZDZ.closest.MCCC))
  #FZDZ.closest.SRPC[ind.FZDZ.closest.MCCC.isNA] <- NA
  # calculate some statistics
  FZDZ.closest.MCCC.mean      <- mean(FZDZ.closest.MCCC, na.rm=TRUE)
  FZDZ.closest.SRPC.mean      <- mean(FZDZ.closest.SRPC, na.rm=TRUE)

  # set SLW.closest.SRPC to NA where AC was not in return
  #ind.SLW.closest.MCCC.isNA  <- which(is.na(SLW.closest.MCCC))
  #SLW.closest.SRPC[ind.SLW.closest.MCCC.isNA] <- NA
  # calculate some statistics
  SLW.closest.MCCC.mean       <- mean(SLW.closest.MCCC, na.rm=TRUE)
  SLW.closest.SRPC.mean       <- mean(SLW.closest.SRPC, na.rm=TRUE)
  print(SLW.closest.MCCC.mean)
  print(SLW.closest.SRPC.mean)

  # set FZDZ.closest.SRPC to NA where AC was not in return
  #ind.MIXPHA.closest.MCCC.isNA<- which(is.na(MIXPHA.closest.MCCC))
  #MIXPHA.closest.SRPC[ind.MIXPHA.closest.MCCC.isNA] <- NA
  # calculate some statistics
  MIXPHA.closest.MCCC.mean    <- mean(MIXPHA.closest.MCCC, na.rm=TRUE)
  MIXPHA.closest.SRPC.mean    <- mean(MIXPHA.closest.SRPC, na.rm=TRUE)

  #------------------------------------
  # output arrays as csv files
  #------------------------------------
  write.csv(FZDZ.closest.MCCC,   file = paste(closest.int.val.array.dir, yyyymmdd.for.case, "_FZDZ.clo.MCCC.csv",   sep=""), row.names=FALSE)
  write.csv(SLW.closest.MCCC,    file = paste(closest.int.val.array.dir, yyyymmdd.for.case, "_SLW.clo.MCCC.csv",    sep=""), row.names=FALSE)
  write.csv(MIXPHA.closest.MCCC, file = paste(closest.int.val.array.dir, yyyymmdd.for.case, "_MIXPHA.clo.MCCC.csv", sep=""), row.names=FALSE)
  write.csv(FZDZ.closest.SRPC,   file = paste(closest.int.val.array.dir, yyyymmdd.for.case, "_FZDZ.clo.SRPC.csv",   sep=""), row.names=FALSE)
  write.csv(SLW.closest.SRPC,    file = paste(closest.int.val.array.dir, yyyymmdd.for.case, "_SLW.clo.SRPC.csv",    sep=""), row.names=FALSE)
  write.csv(MIXPHA.closest.SRPC, file = paste(closest.int.val.array.dir, yyyymmdd.for.case, "_MIXPHA.clo.SRPC.csv", sep=""), row.names=FALSE)
  
  # SETUP READ OF CSV ARRAY FILES HERE, SO CLOSEST SEARCH AND FILE LOAD DOES NOT NEED REPEATING IF NO PREVIOUS CHANGES HAVE BEEN MADE
  
  #------------------------------------
  # plots
  #------------------------------------
  remove(ind.hhmm.closest.wval)
  ind.hhmm.closest.wval        <- which(!is.na(hhmm.closest))
  hhmm.closest.indomain        <- as.POSIXct(flight.track.df$ac.date.posix[ind.hhmm.closest.wval])
  FZDZ.closest.MCCC.indomain   <- FZDZ.closest.MCCC[ind.hhmm.closest.wval]
  FZDZ.closest.SRPC.indomain   <- FZDZ.closest.SRPC[ind.hhmm.closest.wval]
  SLW.closest.MCCC.indomain    <- SLW.closest.MCCC[ind.hhmm.closest.wval]
  SLW.closest.SRPC.indomain    <- SLW.closest.SRPC[ind.hhmm.closest.wval]
  MIXPHA.closest.MCCC.indomain <- MIXPHA.closest.MCCC[ind.hhmm.closest.wval]
  MIXPHA.closest.SRPC.indomain <- MIXPHA.closest.SRPC[ind.hhmm.closest.wval]
  
  dat.indomain.df              <- data.frame(hhmm.closest.indomain, FZDZ.closest.MCCC.indomain, SLW.closest.MCCC.indomain, MIXPHA.closest.MCCC.indomain, FZDZ.closest.SRPC.indomain, SLW.closest.SRPC.indomain, MIXPHA.closest.SRPC.indomain)
  colnames(dat.indomain.df)    <- c("HH:MM", "FZDZ int (mosaic)", "SLW int (mosaic)", "MIXPHA int (mosaic)", "FZDZ int (polar)", "SLW int (polar)", "MIXPHA int (polar")
  
  #FZDZ 
  #  MCCC and SRPC timeseries
  par(mfrow=c(1,2))
  par(mar=c(4, 4, 3, 3))
  image(t(flip(FZDZ.MCCC.raster, direction="x")), xlim=c(min(flight.track.df$lon)-0.5, max(flight.track.df$lon)+0.5), ylim=c(min(flight.track.df$lat)-0.5, max(flight.track.df$lat)+0.5), xlab="Lon [deg]", ylab="Lat [deg]", main="FZDZ INT, MCCC, 1.0km")
  #image(t(flip(FZDZ.MCCC.raster, direction="x")), xlim=c(flight.track.df$lon[min(ind.flttrk.times)]-1.5, flight.track.df$lon[max(ind.flttrk.times)]+1.5), ylim=c(flight.track.df$lat[min(ind.flttrk.times)-1.5], flight.track.df$lat[max(ind.flttrk.times)]+1.5), xlab="Lon [deg]", ylab="Lat [deg]", main="FZDZ INT, MCCC, 1.0km")
  lines(flight.track.df$lon[ind.flttrk.times], flight.track.df$lat[ind.flttrk.times])
  text(radar.lon.deg,     radar.lat.deg,     "*",         cex=2, col="blue")
  text(radar.lon.deg+0.0, radar.lat.deg+0.05, radar.for.case,     col="blue")
  #text(flight.track.df$lon[620], flight.track.df$lat[620], "*", cex=3)
  grid()
  image(t(flip(FZDZ.SRPC.raster, direction="x")), xlim=c(-270,-460), ylim=c(300, 360), xlab="Range [km]", ylab="Azim [deg]", main="FZDZ INT, SRPC, 0.5deg")
  grid()
  ggplot() +
    #geom_point(aes(color = as.factor(category))) + 
    geom_line( data=dat.indomain.df, aes(x = hhmm.closest.indomain, y = FZDZ.closest.MCCC.indomain), color="blue") +
    geom_line( data=dat.indomain.df, aes(x = hhmm.closest.indomain, y = FZDZ.closest.SRPC.indomain), color="red") +
    geom_point(data=dat.indomain.df, aes(x = hhmm.closest.indomain, y = FZDZ.closest.MCCC.indomain), color="blue") +
    geom_point(data=dat.indomain.df, aes(x = hhmm.closest.indomain, y = FZDZ.closest.SRPC.indomain), color="red") +
    ylim(0, 1) +
    labs(x="HH:MM [UTC]", y="FZDZ int") +
    scale_x_datetime(date_labels = "%H:%M")
  #plot(hhmm.closest.indomain, FZDZ.closest.MCCC.indomain, ylim=c(0,1), xlab="Time", ylab="FZDZ INT")
  #text(400, 0.4, paste("FZDZ MCCC mean = ", round(FZDZ.closest.MCCC.mean, digits=2), sep=""))
  #grid()
  #plot(hhmm.closest.indomain, FZDZ.closest.SRPC.indomain, ylim=c(0,1), xlab="Time", ylab="FZDZ INT")
  #text(400, 0.4, paste("FZDZ SRPC mean = ", round(FZDZ.closest.SRPC.mean, digits=2), sep=""))
  #grid()

  # SLW MCCC and SRPC timeseries
  par(mfrow=c(2,2))
  par(mar=c(4, 4, 3, 3))
  image(t(flip(SLW.MCCC.raster, direction="x")), xlim=c(lon.flttrk.mean-1.5, lon.flttrk.mean+1.5), ylim=c(lat.flttrk.mean-1.5, lat.flttrk.mean+1.5), xlab="Lon [deg]", ylab="Lat [deg]", main="SLW INT, MCCC, 1.0km")
  lines(flight.track.df$lon, flight.track.df$lat)
  text(radar.lon.deg,     radar.lat.deg,     "*",         cex=2, col="blue")
  text(radar.lon.deg+0.0, radar.lat.deg+0.1, radar.for.case,     col="blue")
  grid()
  image(t(flip(SLW.SRPC.raster, direction="x")), xlim=c(-300, -458), ylim=c(230, 360), xlab="Range [km]", ylab="Azim [deg]", main="SLW INT, SRPC, 0.5deg")
  grid()
  plot(SLW.closest.MCCC.indomain, ylim=c(0,1))
  grid()
  ggplot() +
    #geom_point(aes(color = as.factor(category))) + 
    geom_line(data=dat.indomain.df,  aes(x = hhmm.closest.indomain, y = SLW.closest.MCCC.indomain), color="blue") +
    geom_point(data=dat.indomain.df, aes(x = hhmm.closest.indomain, y = SLW.closest.MCCC.indomain), color="blue") +
    geom_line(data=dat.indomain.df,  aes(x = hhmm.closest.indomain, y = SLW.closest.SRPC.indomain), color="red") +
    geom_point(data=dat.indomain.df, aes(x = hhmm.closest.indomain, y = SLW.closest.SRPC.indomain), color="red") +
    ylim(0, 1) +
    labs(x="HH:MM [UTC]", y="SLW int") +
    scale_x_datetime(date_labels = "%H:%M")
  #text(400, 0.4, paste("SLW MCCC mean = ", round(SLW.closest.MCCC.mean, digits=2), sep=""))
  #grid()
  #plot(SLW.closest.SRPC.indomain, ylim=c(0,1), xlab="Time", ylab="SLW INT")
  #text(400, 0.4, paste("SLW SRPC mean = ", round(SLW.closest.SRPC.mean, digits=2), sep=""))
  #grid()
  
  # MIXPHA MCCC and SRPC timeseries
  par(mfrow=c(2,2))
  par(mar=c(4, 4, 3, 3))
  image(t(flip(MIXPHA.MCCC.raster, direction="x")), xlim=c(lon.flttrk.mean-1.5, lon.flttrk.mean+1.5), ylim=c(lat.flttrk.mean-1.5, lat.flttrk.mean+1.5), xlab="Lon [deg]", ylab="Lat [deg]", main="MIXPHA INT, MCCC, 1.0km")
  lines(flight.track.df$lon, flight.track.df$lat)
  text(radar.lon.deg,     radar.lat.deg,     "*",         cex=2, col="blue")
  text(radar.lon.deg+0.0, radar.lat.deg+0.1, radar.for.case,     col="blue")
  grid()
  image(t(flip(MIXPHA.SRPC.raster, direction="x")), xlim=c(-300, -458), ylim=c(230, 360), xlab="Range [km]", ylab="Azim [deg]", main="MIXPHA INT, SRPC, 0.5deg")
  grid()
  plot(MIXPHA.closest.MCCC.indomain, ylim=c(0,1))
  grid()
  ggplot() +
    #geom_point(aes(color = as.factor(category))) + 
    geom_line(data=dat.indomain.df,  aes(x = hhmm.closest.indomain, y = MIXPHA.closest.MCCC.indomain), color="blue") +
    geom_point(data=dat.indomain.df, aes(x = hhmm.closest.indomain, y = MIXPHA.closest.MCCC.indomain), color="blue") +
    geom_line(data=dat.indomain.df,  aes(x = hhmm.closest.indomain, y = MIXPHA.closest.SRPC.indomain), color="red") +
    geom_point(data=dat.indomain.df, aes(x = hhmm.closest.indomain, y = MIXPHA.closest.SRPC.indomain), color="red") +
    ylim(0, 1) +
    labs(x="HH:MM [UTC]", y="MIXPHA int") +
    scale_x_datetime(date_labels = "%H:%M")
  #text(400, 0.4, paste("MIXPHA MCCC mean = ", round(MIXPHA.closest.MCCC.mean, digits=2), sep=""))
  #grid()
  #plot(MIXPHA.closest.SRPC.indomain, ylim=c(0,1), xlab="Time", ylab="MIXPHA INT")
  #text(400, 0.4, paste("MIXPHA SRPC mean = ", round(MIXPHA.closest.SRPC.mean, digits=2), sep=""))
  #grid()
  
  
  # FZDZMCCC versus SRPC (scatter)
  ggplot() +
    #geom_point(aes(color = as.factor(category))) + 
    geom_point(data=dat.indomain.df, aes(x = FZDZ.closest.SRPC.indomain, y = FZDZ.closest.MCCC.indomain), color="black") +
    ylim(0, 1) +
    xlim(0, 1) +
    geom_abline(intercept = 0, slope = 1) +
    labs(x="FZDZ int (polar)", y="FZDZ int (mosaic)")
  # SLW MCCC versus SRPC (scatter)
  ggplot() +
    #geom_point(aes(color = as.factor(category))) + 
    geom_point(data=dat.indomain.df, aes(x = SLW.closest.SRPC.indomain, y = SLW.closest.MCCC.indomain), color="black") +
    ylim(0, 1) +
    xlim(0, 1) +
    geom_abline(intercept = 0, slope = 1) +
    labs(x="SLW int (polar)", y="SLW int (mosaic)")
  # MIXPHA MCCC versus SRPC (scatter)
  ggplot() +
    #geom_point(aes(color = as.factor(category))) + 
    geom_point(data=dat.indomain.df, aes(x = MIXPHA.closest.SRPC.indomain, y = MIXPHA.closest.MCCC.indomain), color="black") +
    ylim(0, 1) +
    xlim(0, 1) +
    geom_abline(intercept = 0, slope = 1) +
    labs(x="MIXPHA int ((polar)", y="MIXPHA int (mosaic)")

} #  end of for (i in 1:dim(icicle.df)[1])

## basic stats for FZDZ
#mean(FZDZ.MCCC[260:460, 500:650, 3],     na.rm=TRUE)
#median(FZDZ.MCCC[260:460, 500:650, 3],   na.rm=TRUE)
#mean(FZDZ.SRPC[1:733, 432:720, 1],       na.rm=TRUE)
#median(FZDZ.SRPC[1:733, 432:720, 1],     na.rm=TRUE)
## basic stats for SLW
#mean(SLW.MCCC[260:460, 500:650, 3],      na.rm=TRUE)
#median(SLW.MCCC[260:460, 500:650, 3],    na.rm=TRUE)
#mean(SLW.SRPC[1:733, 432:720, 1],        na.rm=TRUE)
#median(SLW.SRPC[1:733, 432:720, 1],      na.rm=TRUE)
## basic stats for MIXPHA
#mean(MIXPHA.MCCC[260:460, 500:650, 3],   na.rm=TRUE)
#median(MIXPHA.MCCC[260:460, 500:650, 3], na.rm=TRUE)
#mean(MIXPHA.SRPC[1:733, 432:655, 1],     na.rm=TRUE)
#median(MIXPHA.SRPC[1:733, 432:655, 1],   na.rm=TRUE)

# OLDER STUFF BELOW.  ABOVE CODED 11/5/2019
##----------------------------
## build lat/lon arrays and matrices
##----------------------------
## for SRPC
#len.FZDZ.SRPC.range        <- dim(FZDZ.SRPC)[1]
#len.FZDZ.SRPC.azim         <- dim(FZDZ.SRPC)[2]
#len.FZDZ.SRPC.tilt         <- dim(FZDZ.SRPC)[3]
#RadIA.SRPC.range.atom.km   <- seq(from = 2.125, to = 459.875, by = 0.25)                              # length is 1832
#RadIA.SRPC.range.matr.km   <- rep.col(RadIA.SRPC.range.atom.km, 720)
#RadIA.SRPC.azim.atom.deg   <- seq(from = 0.0, to = 359.5,  by = (359.5-0.0)/(len.FZDZ.SRPC.azim - 1)) # length is 720
#RadIA.SRPC.azim.matr.deg   <- rep.row(RadIA.SRPC.azim.atom.deg, 1832)
##RadIA.SRPC.tilt            <- c()

## for MCCC
#len.FZDZ.MCCC              <- dim(FZDZ.MCCC)[1]
#RadIA.MCCC.lon             <- seq(from = ICICLE.max.lon, to = ICICLE.min.lon, by = -(ICICLE.max.lon-ICICLE.min.lon)/(len.FZDZ.MCCC - 1))
#RadIA.MCCC.lat             <- seq(from = ICICLE.min.lat, to = ICICLE.max.lat, by =  (ICICLE.max.lat-ICICLE.min.lat)/(len.FZDZ.MCCC - 1))

##----------------------------
## calculate alt of beam at each tilt and range, based on tilt angle and radar alt
##----------------------------
## convert NEXRAD volume tilt angles from deg to rad
#theta.rad                  <- NISTdegTOradian(theta.deg)

#num.theta.deg              <- length(theta.rad)

## make a matrix of arrays, one for each tilt angle
#alt.beam.km                <- replicate(num.theta.deg, RadIA.SRPC.range.matr.km)   # make a matrix for every tilt
#for (j in 1:length(theta.rad)) {
#  alt.beam.km[,,j]         <- (RadIA.SRPC.range.matr.km * sin(theta.rad[j])) + ((RadIA.SRPC.range.matr.km^2)/ radius.earth.km)
#}

##----------------------------
## print some output to GUI
##----------------------------
#print(paste("MCCC lon value for profile: ", round(RadIA.MCCC.lon[ind.lon], digits=3)), sep="")
#print(paste("MCCC lat value: for profile", round(RadIA.MCCC.lat[ind.lat], digits=3)), sep="")

##----------------------------
## plots
##----------------------------
## FZDZ.SRPC
#col.l                      <- colorRampPalette(c('purple', 'blue', 'green', 'yellow', 'red'))(30)
#levelplot(FZDZ.SRPC[ , , ind.SRPC.alt], xlim=c(1, 600), ylim=c(400, 720), col.regions=col.l, xlab="range", ylab="azim")

## FZDZ.MCCC
#col.l                      <- colorRampPalette(c('purple', 'blue', 'green', 'yellow', 'red'))(30)
#levelplot(FZDZ.MCCC[ , ,ind.nssl.mosaic.alt], nssl.mosaic.alt.array, xlim=c(250, 500), ylim=c(450, 650), col.regions=col.l, xlab="lon [deg]", ylab="lat [deg]", cuts=10, scales=list(x=list(at=c(250,375,500), labels=c(round(RadIA.MCCC.lon[500], digits=2), round(RadIA.MCCC.lon[375], digits=2), round(RadIA.MCCC.lon[250], digits=2))), y=list(at=c(450,550,650), labels=c(round(RadIA.MCCC.lat[450], digits=2), round(RadIA.MCCC.lat[550], digits=2), round(RadIA.MCCC.lat[650], digits=2)))), main=paste(radia.date, " ", substring(radia.MCCC.filename, 10,11), ":", substring(radia.MCCC.filename, 12,13), "Z, ", nssl.mosaic.alt.array[ind.nssl.mosaic.alt],"km, RadIA FZDZ INT", sep="")) + layer(panel.text(360, 600, "x"))

## plot profile of given RadIA MCCC algo interest at a given location
#par(mfrow=c(1,3))
#par(mar=c(7, 4, 2, 1.5))
#plot(FZDZ.MCCC[360, 600, ], nssl.mosaic.alt.array, pch=20, type="b", xlim=c(0, 1), col="blue", lwd=2, ylab="Alt [km MSL]", main=paste(radia.date, " FZDZ int", sep=""))
#abline(v=0.40, lty=2, lwd=2)
#grid()
#plot(SLW.MCCC[360, 600, ], nssl.mosaic.alt.array, pch=20, type="b", xlim=c(0, 1), col="green", lwd=2, ylab="Alt [km MSL]", main=paste(radia.date, " SLW int", sep=""))
#abline(v=0.76, lty=2, lwd=2)
#grid()
#plot(MIXPHA.MCCC[360, 600, ], nssl.mosaic.alt.array, pch=20, type="b", xlim=c(0, 1), col="orange", lwd=2, ylab="Alt [km MSL]", main=paste(radia.date, " MIXPHA int", sep=""))
#abline(v=0.42, lty=2, lwd=2)
#grid()

#par(mfrow=c(1,1))
#hh.mm.str.F24              <- substr(as.character(flight.track.df$ac.date.posix), 12, 16)
#mean.ICICLE.flight.lat.F24 <- mean(flight.track.df$lat)
#mean.ICICLE.flight.lon.F24 <- mean(flight.track.df$lon)
#ICICLE.F24.map             <- get_map(location = c(mean.ICICLE.flight.lon.F24-1.0, mean.ICICLE.flight.lat.F24+1.0), maptype = "terrain", source = "google", zoom = 7)
#ggmap(ICICLE.F24.map)
##   create range rings from each user defined radar
#ind.radar                  <- which(NEXRAD.site.df$ICAO == " KDVN")
#datatest.df                <- data.frame(ID = as.numeric(c(1:3)), longitude = as.numeric(c(NEXRAD.site.df$lon[38], NEXRAD.site.df$lon[35], NEXRAD.site.df$lon[7])), latitude = as.numeric(c(NEXRAD.site.df$lat[38], NEXRAD.site.df$lat[35], NEXRAD.site.df$lat[7])))
#myCircles.125km            <- make_circles(datatest.df, 125.0)
#myCircles.100km            <- make_circles(datatest.df, 100.0)
#myCircles.050km            <- make_circles(datatest.df, 50.0)
#ICICLE.F24.map             <- get_map(location = c(mean.ICICLE.flight.lon.F24-1.0, mean.ICICLE.flight.lat.F24+1.0), maptype = "terrain", source = "google", zoom = 7)
#ggmap(ICICLE.F24.map)
#ggmap(ICICLE.F24.map) + 
#  geom_point(data = flight.track.df, aes(x = lon, y = lat), alpha = 0.5, size = 0.5, color="green") +
#  #geom_point(data = ICICLEflighttrack.F24, aes(x = lon, y = lat), alpha = 0.5, size = 0.5, color="black") +
#  labs(title = "ICICLE flight track for F24 2019-02-26") +
#  #annotate("text", x=ICICLEflighttrack.F24$lon[150]+0.00, y=ICICLEflighttrack.F24$lat[150]+0.05,label=paste(hh.mm.str.F24[150], "<", sep=""),      color="blue") +
#  #annotate("text", x=ICICLEflighttrack.F24$lon[338]+0.00, y=ICICLEflighttrack.F24$lat[338]-0.05,label=paste(hh.mm.str.F24[338], "FZDZ<", sep=""),      color="blue") +
#  #annotate("text", x=ICICLEflighttrack.F24$lon[389]+0.00, y=ICICLEflighttrack.F24$lat[389]-0.10,label=paste(hh.mm.str.F24[389],  "MA", sep=""),        color="blue") +
#  #annotate("text", x=ICICLEflighttrack.F24$lon[516]+0.00, y=ICICLEflighttrack.F24$lat[516]-0.05,label=paste(hh.mm.str.F24[516],  "", sep=""),          color="blue") +
#  #annotate("text", x=ICICLEflighttrack.F24$lon[531]+0.00, y=ICICLEflighttrack.F24$lat[531]+0.08,label=paste(hh.mm.str.F24[531],  "FZDZ/SLD<", sep=""), color="blue") +
#  #annotate("text", x=ICICLEflighttrack.F24$lon[581]+0.00, y=ICICLEflighttrack.F24$lat[581]+0.03,label=paste(hh.mm.str.F24[581],  "SLD", sep=""),       color="blue") +
#  #annotate("text", x=ICICLEflighttrack.F24$lon[587]+0.00, y=ICICLEflighttrack.F24$lat[587]+0.04,label=paste(hh.mm.str.F24[587],  "FZDZ", sep=""),       color="blue") +
#  #annotate("text", x=ICICLEflighttrack.F24$lon[623]+0.00, y=ICICLEflighttrack.F24$lat[623]-0.10,label=paste(hh.mm.str.F24[623],  "FZDZ", sep=""),      color="blue") +
#  #annotate("text", x=ICICLEflighttrack.F24$lon[662]+0.00, y=ICICLEflighttrack.F24$lat[662]-0.05,label=paste(hh.mm.str.F24[662],  "FZDZ>", sep=""),     color="blue") +
#  annotate("text", x=NEXRAD.site.df$lon[38],              y=NEXRAD.site.df$lat[38],             label = "*", size = 3,                                 color="red") +
#  annotate("text", x=NEXRAD.site.df$lon[38]+0.10,         y=NEXRAD.site.df$lat[38]+0.10,        label = NEXRAD.site.df$ICAO[38], size = 4,             color="red") +
#  annotate("text", x=NEXRAD.site.df$lon[35],              y=NEXRAD.site.df$lat[35],             label = "*", size = 3,                                 color="red") +
#  annotate("text", x=NEXRAD.site.df$lon[35]+0.10,         y=NEXRAD.site.df$lat[35]+0.10,        label = NEXRAD.site.df$ICAO[35], size = 4,             color="red") +
#  annotate("text", x=NEXRAD.site.df$lon[7],               y=NEXRAD.site.df$lat[7],              label = "*", size = 3,                                 color="red") +
#  annotate("text", x=NEXRAD.site.df$lon[7]+0.10,          y=NEXRAD.site.df$lat[7]+0.10,         label = NEXRAD.site.df$ICAO[7], size = 4,              color="red") +
#  geom_polygon(data = myCircles.125km, aes(lon, lat, group = ID), color = "black", alpha = 0) +
#  geom_polygon(data = myCircles.100km, aes(lon, lat, group = ID), color = "black", alpha = 0) +
#  geom_polygon(data = myCircles.050km, aes(lon, lat, group = ID), color = "black", alpha = 0)
##grid()

## convert FZDZ.MCCC to raster
#dim(FZDZ.MCCC.df)
#RadIA.lat.df        <- data.frame(seq(from=36, to=48, by=(48-36)/(dim(FZDZ.MCCC)[1]-1)))
#RadIA.lon.df        <- data.frame(seq(from=83, to=95, by=(95-83)/(dim(FZDZ.MCCC)[1]-1)))
#names(RadIA.lat.df) <- "lat"
#names(RadIA.lon.df) <- "lon"
#FZDZ.MCCC.new.df    <- rbind(FZDZ.MCCC.df, RadIA.lat.df)

##
#usa         <- map_data("usa")
#states      <- map_data("state")
#ICICLE.FL24 <- subset(states, region %in% c("iowa"))
#counties    <- map_data("county")
#ia_county   <- subset(counties, region == "iowa")
##ggplot(data=states) + geom_polygon(aes(x=long, y = lat, fill=region, group = group), color="white") + 
##  coord_fixed(1.3) +
##  guides(fill=FALSE)

##
#ia_base     <- ggplot(data = ICICLE.FL24, mapping = aes(x = long, y = lat, group = group)) + 
#  coord_fixed(1.3) + 
#  geom_polygon(color = "black", fill = "gray")
#ia_base + theme_nothing() + 
#  geom_polygon(data = ia_county, fill = NA, color = "white") +
#  geom_polygon(color = "black", fill = NA)  +
#  annotate("text", x=NEXRAD.site.df$lon[38],              y=NEXRAD.site.df$lat[38],             label = "*", size = 3,                                 color="red") +
#  annotate("text", x=NEXRAD.site.df$lon[38]+0.10,         y=NEXRAD.site.df$lat[38]+0.10,        label = NEXRAD.site.df$ICAO[38], size = 4,             color="red") +
#  annotate("text", x=NEXRAD.site.df$lon[35],              y=NEXRAD.site.df$lat[35],             label = "*", size = 3,                                 color="red") +
#  annotate("text", x=NEXRAD.site.df$lon[35]+0.10,         y=NEXRAD.site.df$lat[35]+0.10,        label = NEXRAD.site.df$ICAO[35], size = 4,             color="red") +
#  annotate("text", x=NEXRAD.site.df$lon[7],               y=NEXRAD.site.df$lat[7],              label = "*", size = 3,                                 color="red") +
#  annotate("text", x=NEXRAD.site.df$lon[7]+0.10,          y=NEXRAD.site.df$lat[7]+0.10,         label = NEXRAD.site.df$ICAO[7], size = 4,              color="red") +
#  geom_polygon(data = myCircles.125km, aes(lon, lat, group = ID), color = "black", alpha = 0) +
#  geom_polygon(data = myCircles.100km, aes(lon, lat, group = ID), color = "black", alpha = 0) +
#  geom_polygon(data = myCircles.050km, aes(lon, lat, group = ID), color = "black", alpha = 0)

##
#data(wrld_simpl)
#data(usaMapEnv)
#r      <- raster(wrld_simpl, res=1)
#wrld_r <- rasterize(wrld_simpl, r)
#par(mfrow=c(1,1))
#plot(wrld_r, col = "grey")
#plot(raster(FZDZ.MCCC[,,3]), add = TRUE)

#ICICLE.FL24.map <- raster(ICICLE.FL24)

## alaska example from: http://zevross.com/blog/2015/03/30/map-and-analyze-raster-data-in-r/
#url         <- "http://qgis.org/downloads/data/qgis_sample_data.zip"
#mydir       <- "/d1/serke/projects/RADIA_FAA"
#temp        <- tempfile(tmpdir=mydir, fileext=".zip")
#download.file(url, temp)
#unzip(temp, exdir=mydir)
#unlink(temp) #delete the zip file
#fpath       <- list.files(path = mydir, full.names = TRUE, pattern = "qgis_sample_data")
##fpath       <- gsub("/", "\\\\", fpath)
#fpath
#landusepath <- paste(fpath, "/raster/landcover.img", sep="")
#landusepath
#landuse.raw <- raster(landusepath)
#plot(landuse.raw, axes=FALSE)
#vals        <- unique(values(landuse.raw))
#recl        <- matrix(c(vals, c(0, rep(1, 6), 9, 1,1, 13)),ncol=2)
#recl
#landuse     <- reclassify(landuse.raw, rcl=recl)
#plot(landuse, legend=FALSE, axes=FALSE)
## Regions polygon shapefile
#regionpath  <- paste(fpath, "/shapefiles", sep="")
#region      <- readOGR(dsn=regionpath, layer="regions") 
## we will use ggplot to plot the regions
#ggplot() + geom_polygon(data=region,  aes(x=long, y=lat, group=group), fill="cadetblue", color="grey") +
#           coord_equal()+xlim(c(-5000000, 5000000))+ylim(c(1000000, 8000000))
## Create a subset with our regions of interest
#myregions   <- c( "Anchorage", "Yukon-Koyukuk", "North Slope")
#region.sm   <- region[region$NAME_2 %in% myregions,]
## crop, rasterize and mask 
#cr          <- crop(landuse, region.sm)
#fr          <- rasterize(region.sm, cr)
#lr          <- mask(x=cr, mask=fr)
## let's map those pieces so you can see the result. Since I just
## want the raster with no legend/axes etc I'm creating a function
## to strip the plot
#nakedMap    <- function(dat, title="") {
#  gplot(dat)+geom_tile(aes(fill=value))+
#    ggtitle(title)+
#    coord_equal()+ 
#    scale_x_continuous(expand = c(0,0)) + 
#    scale_y_continuous(expand = c(0,0)) +
#    theme(line = element_blank(),
#          line = element_blank(),
#          axis.text=element_blank(),
#          axis.title=element_blank(),
#          legend.position = "none")
#}
#cr.plot     <- nakedMap(cr, title="Cropped")
#fr.plot     <- nakedMap(fr, title="Regions")
#lr.plot     <- nakedMap(lr, title="Masked")
#grid.arrange(cr.plot, fr.plot, lr.plot, ncol=3) #use package gridExtra
## centroids for the labels
#centroids   <- cbind(coordinates(region.sm), region.sm@data)
#names(centroids)[1:2]<-c("x", "y")
## use gplot (not ggplot) from rasterVis
## geom_tile adds the raster, geom_polygon adds the regions
## geom_text adds the labels at the centroids
#gplot(lr)+
#  geom_tile(aes(fill=factor(value, labels=c("Water", "Green", "Shrubland", "Urban"))), alpha=0.8)+
#  scale_fill_manual(values = c("steelblue3", "forestgreen", "ghostwhite", "red"),
#                    name= "Land use code")+
#  geom_polygon(data=region.sm, aes(x=long, y=lat, group=group), 
#               fill=NA,color="grey50", size=1)+
#  geom_text(data=centroids, aes(x=x, y=y, label=NAME_2), fontface="bold")+
#  coord_equal()
## Elevation raster
#elevpath    <- paste(fpath, "/raster/SR_50M_alaska_nad.tif", sep="")
#elev        <- raster(elevpath)
## Are the resolutions the same
#res(landuse.raw)
### [1] 3280 3280
#res(elev)
### [1] 7181.354 7181.354
## No, we will use resample to make the same resolution
## NOTE: this takes a minute or so to run
#elev        <- resample(elev, landuse.raw,method='ngb')
## Check again to see that the resolutions match
#res(landuse.raw)
### [1] 3280 3280
#res(elev)
### [1] 3280 3280

## we manually selected an extent to crop to
#ext<-extent(-2500000, 2500000, 2000000, 7650000)
#elev<-crop(elev, ext)
#landuse.raw<-crop(landuse.raw, ext)