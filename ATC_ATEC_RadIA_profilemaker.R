#----------------------------------------------
#
# Name:    ATC_ATEC_RadIA_profilemaker.R
#
# Purpose: output profiles of RadIA interests over ATC, MD for further product development
#
# Created: 4.25.2019 dserke
#
#----------------------------------------------

# add libraries
library(lubridate)
library(R.matlab)
library(ncdf4)
library(NISTunits)
library(useful)

require(dplyr)
require(geosphere)

# remove scientific notation in printing
options(scipen=999)

#----------------------------------------------
# constants and parameters
#----------------------------------------------
# define yyyymmdd character string to work from
radia.yyyy                         <- "2019"
radia.mm                           <- "02"
radia.dd                           <- "12"
radia.date                         <- paste(radia.yyyy, radia.mm, radia.dd, sep = "")

radius.earth.km                    <- 6378.1

tilts.deg                          <- c(0.48, 0.87, 1.32, 1.80, 2.42, 3.12, 4.00, 5.10, 6.42)
tilts.rad                          <- NISTdegTOradian(tilts.deg)

# distance from KDOX to ATC (measured from Google Maps)
dist.min.km        <- 87.5
dist.max.km        <- 96.5

# azimuth from KDOX to ATC (measured from CIDD display)
azim.min.deg       <- 308
azim.max.deg       <- 323
ind.azim.min       <- 2 * azim.min.deg
ind.azim.max       <- 2 * azim.max.deg

#----------------------------------------------
# define data paths
#----------------------------------------------
# define RadIA SRPC nc interest directories
radia.SRPC.nc.FZDZ.dir             <- file.path(paste("/d1/serke/projects/case_studies/ATEC_ranges_2019/Aberdeen_MD/data/RadIA/nc/FZDZ/",   radia.date, "/", sep=""))
radia.SRPC.nc.MIXPHA.dir           <- file.path(paste("/d1/serke/projects/case_studies/ATEC_ranges_2019/Aberdeen_MD/data/RadIA/nc/MIXPHA/", radia.date, "/", sep=""))
radia.SRPC.nc.SLW.dir              <- file.path(paste("/d1/serke/projects/case_studies/ATEC_ranges_2019/Aberdeen_MD/data/RadIA/nc/SLW/",    radia.date, "/", sep=""))
#radia.SRPC.nc.PLATES.dir           <- file.path(paste("/d1/serke/projects/case_studies/ATEC_ranges_2019/Aberdeen_MD/data/RadIA/nc/PLATES/", radia.date, "/", sep=""))

# define NEXRAD site location csv file path
nexrad.site.dataframetxt.dir       <- file.path("/d1/serke/projects/case_studies/SNOWIE/data/RadIA_data/nexrad_site_data/")

# output data
#output.file.name                   <- file.path("/d1/serke/projects/case_studies/SNOWIE/data/RadIA_matchto_AC_data/out.csv")

#----------------------------------------------
# define a list of nc files in dir
#----------------------------------------------
# for FZDZ
radia.SRPC.nc.FZDZ.filelist.str   <- list.files(radia.SRPC.nc.FZDZ.dir,  pattern = "[.nc]",       full.names = FALSE, ignore.case = TRUE)
radia.SRPC.nc.FZDZ.hhmmsslist.str <- substring(radia.SRPC.nc.FZDZ.filelist.str, 1, 6)
radia.SRPC.nc.FZDZ.hhmmsslist.num <- unlist(lapply(radia.SRPC.nc.FZDZ.hhmmsslist.str, as.numeric))
# print results to screen for debugging
radia.SRPC.nc.FZDZ.filelist.str 

# for SLW
radia.SRPC.nc.SLW.filelist.str   <- list.files(radia.SRPC.nc.SLW.dir,  pattern = "[.nc]",       full.names = FALSE, ignore.case = TRUE)
radia.SRPC.nc.SLW.hhmmsslist.str <- substring(radia.SRPC.nc.SLW.filelist.str, 1, 6)
radia.SRPC.nc.SLW.hhmmsslist.num <- unlist(lapply(radia.SRPC.nc.SLW.hhmmsslist.str, as.numeric))

# for MIXPHA
radia.SRPC.nc.MIXPHA.filelist.str   <- list.files(radia.SRPC.nc.MIXPHA.dir,  pattern = "[.nc]",       full.names = FALSE, ignore.case = TRUE)
radia.SRPC.nc.MIXPHA.hhmmsslist.str <- substring(radia.SRPC.nc.MIXPHA.filelist.str, 1, 6)
radia.SRPC.nc.MIXPHA.hhmmsslist.num <- unlist(lapply(radia.SRPC.nc.MIXPHA.hhmmsslist.str, as.numeric))

#----------------------------------------------
# load the NEXRAD site location text file
#----------------------------------------------
NEXRAD.site.df                    <- read.csv(paste(nexrad.site.dataframetxt.dir, "nexrad_site.csv", sep = ""), header = FALSE, sep = ",", dec = ".", stringsAsFactors=FALSE)
colnames(NEXRAD.site.df)          <- c("NCDCID", "ICAO", "WBAN", "radname", "COUNTRY", "STATE", "COUNTY", "lat", "lon", "elev", "GMTdiff", "STN_TYPE")
head(NEXRAD.site.df)

# load info specific to KDOX
ind.KDOX                          <- which(NEXRAD.site.df$ICAO == " KDOX")
KDOX.elev.ft                      <- NEXRAD.site.df$elev[ind.KDOX]
KDOX.lat.deg                      <- NEXRAD.site.df$lat[ind.KDOX]
KDOX.lon.deg                      <- NEXRAD.site.df$lon[ind.KDOX]

#----------------------------------------------
# loop through all of the available nc files
#----------------------------------------------
for (k in 1:length(radia.SRPC.nc.FZDZ.filelist.str)) {

  radia.hh <- substring(radia.SRPC.nc.FZDZ.hhmmsslist.str[k], 1, 2)
  radia.mm <- substring(radia.SRPC.nc.FZDZ.hhmmsslist.str[k], 3, 4)
  
  #----------------------------------------------
  # load RadIA interest files 
  #----------------------------------------------
  # load FZDZ nc file
  nc.radia.frzdrz.filename           <- paste(radia.SRPC.nc.FZDZ.dir, radia.SRPC.nc.FZDZ.filelist.str[k], sep="")
  print(paste("loading ", nc.radia.frzdrz.filename, sep=""))
  if (substr(nc.radia.frzdrz.filename, nchar(nc.radia.frzdrz.filename)-1, nchar(nc.radia.frzdrz.filename)) == "gz") {
    print("File is gzipped. Unzipping...")
    untar(nc.radia.frzdrz.filename)
    nc.radia.frzdrz.filename <- substr(nc.radia.frzdrz.filename, 1, nchar(nc.radia.frzdrz.filename)-3)
  }
  nc.radia.frzdrz                    <- nc_open(nc.radia.frzdrz.filename, write = FALSE, verbose = TRUE)
  print(paste("The file has", nc.radia.frzdrz$nvars, "variables"))
  radia.frzdrz.var.num               <- c(1)
  for (i in 1:length(radia.frzdrz.var.num)) {
    radia.frzdrz.nam <- paste("v", radia.frzdrz.var.num[i], sep = "")
    assign(radia.frzdrz.nam, nc.radia.frzdrz$var[[radia.frzdrz.var.num[i]]])
  }
  FZDZ                            <- ncvar_get( nc.radia.frzdrz, v1 )
  print(paste("V1 has name", v1$name))
  nc_close(nc.radia.frzdrz)
  
  # load SLW nc file
  nc.radia.slw.filename           <- paste(radia.SRPC.nc.SLW.dir, radia.SRPC.nc.SLW.filelist.str[k], sep="")
  print(paste("loading ", nc.radia.slw.filename, sep=""))
  if (substr(nc.radia.slw.filename, nchar(nc.radia.slw.filename)-1, nchar(nc.radia.slw.filename)) == "gz") {
    print("File is gzipped. Unzipping...")
    untar(nc.radia.slw.filename)
    nc.radia.slw.filename <- substr(nc.radia.slw.filename, 1, nchar(nc.radia.slw.filename)-3)
  }
  nc.radia.slw                    <- nc_open(nc.radia.slw.filename, write = FALSE, verbose = TRUE)
  print(paste("The file has", nc.radia.slw$nvars, "variables"))
  radia.slw.var.num               <- c(1)
  for (i in 1:length(radia.slw.var.num)) {
    radia.slw.nam <- paste("v", radia.slw.var.num[i], sep = "")
    assign(radia.slw.nam, nc.radia.slw$var[[radia.slw.var.num[i]]])
  }
  SLW                             <- ncvar_get( nc.radia.slw, v1 )
  print(paste("V1 has name", v1$name))
  nc_close(nc.radia.slw)
  
  # load MIXPHA nc file
  nc.radia.mixpha.filename           <- paste(radia.SRPC.nc.MIXPHA.dir, radia.SRPC.nc.MIXPHA.filelist.str[k], sep="")
  print(paste("loading ", nc.radia.mixpha.filename, sep=""))
  if (substr(nc.radia.mixpha.filename, nchar(nc.radia.mixpha.filename)-1, nchar(nc.radia.mixpha.filename)) == "gz") {
    print("File is gzipped. Unzipping...")
    untar(nc.radia.mixpha.filename)
    nc.radia.mixpha.filename <- substr(nc.radia.mixpha.filename, 1, nchar(nc.radia.mixpha.filename)-3)
  }
  nc.radia.mixpha                    <- nc_open(nc.radia.mixpha.filename, write = FALSE, verbose = TRUE)
  print(paste("The file has", nc.radia.mixpha$nvars, "variables"))
  radia.mixpha.var.num               <- c(1)
  for (i in 1:length(radia.mixpha.var.num)) {
    radia.mixpha.nam <- paste("v", radia.mixpha.var.num[i], sep = "")
    assign(radia.mixpha.nam, nc.radia.mixpha$var[[radia.mixpha.var.num[i]]])
  }
  MIXPHA                             <- ncvar_get( nc.radia.mixpha, v1 )
  print(paste("V1 has name", v1$name))
  nc_close(nc.radia.mixpha)
  
  #str(FZDZ)
  
  #----------------------------------------------
  # loop through all of the available volume tilt angles
  #----------------------------------------------
  ATC.FZDZ.max   <- rep(0, dim(FZDZ)[3])
  ATC.FZDZ.med   <- rep(0, dim(FZDZ)[3])
  ATC.FZDZ.ave   <- rep(0, dim(FZDZ)[3])
  ATC.SLW.max    <- rep(0, dim(SLW)[3])
  ATC.SLW.med    <- rep(0, dim(SLW)[3])
  ATC.SLW.ave    <- rep(0, dim(SLW)[3])
  ATC.MIXPHA.max <- rep(0, dim(MIXPHA)[3])
  ATC.MIXPHA.med <- rep(0, dim(MIXPHA)[3])
  ATC.MIXPHA.ave <- rep(0, dim(MIXPHA)[3])
  
  for (j in 1:dim(FZDZ)[3]) {
  
    # convert dist to slant range from KDOX to ATC
    array.range.min      <- dist.min.km / cos(tilts.rad[j])  
    array.range.max      <- dist.max.km / cos(tilts.rad[j]) 
    ind.range.min        <- round((array.range.min - 0.5) / 0.250)
    ind.range.max        <- round((array.range.max - 0.5) / 0.250)
    
    # cut array down to ATC size
    FZDZ.overATC.tiltj   <- FZDZ[ind.range.min:ind.range.max, ind.azim.min:ind.azim.max, j]
    SLW.overATC.tiltj    <- SLW[ind.range.min:ind.range.max, ind.azim.min:ind.azim.max, j]
    MIXPHA.overATC.tiltj <- MIXPHA[ind.range.min:ind.range.max, ind.azim.min:ind.azim.max, j]
    
    # get stats over cut-down array over ATC for FZDZ
    ind.FZDZ.gt0       <- which(FZDZ.overATC.tiltj > 0)
    ATC.FZDZ.max[j]    <- max(FZDZ.overATC.tiltj, na.rm=TRUE)
    if (ATC.FZDZ.max[j] == -Inf) {
      ATC.FZDZ.max[j] <- NA
    }
    ATC.FZDZ.med[j]    <- median(FZDZ.overATC.tiltj[ind.FZDZ.gt0])
    ATC.FZDZ.ave[j]    <- mean(FZDZ.overATC.tiltj[ind.FZDZ.gt0])
    if (is.nan(ATC.FZDZ.ave[j])) {
      ATC.FZDZ.ave[j] <- NA
    }
  
    # get stats over cut-down array over ATC for SLW
    ind.SLW.gt0       <- which(SLW.overATC.tiltj > 0)
    ATC.SLW.max[j]    <- max(SLW.overATC.tiltj, na.rm=TRUE)
    if (ATC.SLW.max[j] == -Inf) {
      ATC.SLW.max[j] <- NA
    }
    ATC.SLW.med[j]    <- median(SLW.overATC.tiltj[ind.SLW.gt0])
    ATC.SLW.ave[j]    <- mean(SLW.overATC.tiltj[ind.SLW.gt0])
    if (is.nan(ATC.SLW.ave[j])) {
      ATC.SLW.ave[j] <- NA
    }

    # get stats over cut-down array over ATC for MIXPHA
    ind.MIXPHA.gt0       <- which(MIXPHA.overATC.tiltj > 0)
    ATC.MIXPHA.max[j]    <- max(MIXPHA.overATC.tiltj, na.rm=TRUE)
    if (ATC.MIXPHA.max[j] == -Inf) {
      ATC.MIXPHA.max[j] <- NA
    }
    ATC.MIXPHA.med[j]    <- median(MIXPHA.overATC.tiltj[ind.MIXPHA.gt0])
    ATC.MIXPHA.ave[j]    <- mean(MIXPHA.overATC.tiltj[ind.MIXPHA.gt0])
    if (is.nan(ATC.MIXPHA.ave[j])) {
      ATC.MIXPHA.ave[j] <- NA
    }
   
    # print some results to the screen
    print("--------------------------------")
    print(paste("k = ", as.character(k), ": Filename = ", radia.SRPC.nc.FZDZ.filelist.str[k]))
    print("--------------------------------")
    print("FZDZ:")
    print(ATC.FZDZ.max)
    print(ATC.FZDZ.med)
    print(ATC.FZDZ.ave)
    print("SLW:")
    print(ATC.SLW.max)
    print(ATC.SLW.med)
    print(ATC.SLW.ave)
    print("MIXPHA:")
    print(ATC.MIXPHA.max)
    print(ATC.MIXPHA.med)
    print(ATC.MIXPHA.ave)
    
    #browser()
    
    # NEED TO MAKE DF HERE:
    # yyyy mm dd hh mm FZDZmax FZDZmed FZDZave SLWmax ......
    RadIA.ints.atATC.df <- data.frame(radia.yyyy, radia.mm, radia.dd, radia.hh, radia.mm, tilts.deg, ATC.FZDZ.max, ATC.FZDZ.med, ATC.FZDZ.ave, ATC.SLW.max, ATC.SLW.med, ATC.SLW.ave, ATC.MIXPHA.max, ATC.MIXPHA.med, ATC.MIXPHA.ave)
    #rm(c(ATC.FZDZ.ave, ATC.FZDZ.max, ATC.FZDZ.med, ATC.MIXPHA.ave, ATC.MIXPHA.max, ATC.MIXPHA.med, ATC.SLW.ave, ATC.SLW.max, ATC.SLW.med))
    
    # AND THEN OUTPUT DF AS CSV
    
    # AND THEN CREATE PROFILE PLOTS OF INTS AND FINAL PRODUCT 
    
  } # end of for (j in 1:length()) {
  
}   # end of for (k in 1:length(radia.SRPC.nc.FZDZ.filelist.str)) {
