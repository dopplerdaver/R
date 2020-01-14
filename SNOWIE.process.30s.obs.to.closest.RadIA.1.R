#----------------------------------------------
#
# Name:    SNOWIE.process.30s.obs.to.closest.RadIA.R
#
# Purpose: 1) load the 30-sec average UWKA flight data plus round particle mat file
#             into a data frame
#          2) load IOP 2ds sizedist ice&liq csv files
#          3) load IOP radia INT files
#          4) output matched 
#
# Created: 4.3.2019 dserke
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
radia.yyyy                         <- "2017"
radia.mm                           <- "01"
radia.dd                           <- "09"
radia.date                         <- paste(radia.yyyy, radia.mm, radia.dd, sep = "")

radius.earth.km                    <- 6378.1

# start and end indices by date for the main loop below

# for 20170108
#start.index <- 78 # for 02:41 to 03:13 Z
#end.index   <- 135
#start.index <- 162 # for 03:26 to 03:13 Z
#end.index   <- 341
#start.index <- 348 # for 05:00 to 05:40 Z
#end.index   <- 429

# for 20170109
#start.index <- 69 # for 04:30 to 04:40 Z
#end.index   <- 89
#start.index <- 89 # for 04:40 to 05:35 Z
#end.index   <- 200
start.index <- 308    # for 06:30 to 07:00 Z
end.index   <- 367

# for 20170118
#start.index <- 80 # for 20:23 to 22:48 Z
#end.index   <- 280
#start.index <- 281 # for 22:48to 22:48 Z
#end.index   <- 370

# for 20170119
#start.index <- 190 # for 16:50 to 17:50 Z
#end.index   <- 310
#start.index <- 204
#end.index   <- 310

# NOTE 174053.nc MIXPHA FILE IS MISSING. CRASHES ON INDEX NUM 286  

# for 20170121
#start.index <- 110 # for 22:40 to 23:00 Z
#end.index   <- 140
#start.index <- 190 # for 23:24 to 23:40 Z
#end.index   <- 220
# NOTE 233720.nc MIXPHA FILE IS MISSING. CRASHES ON INDEX NUM 213 

# for 20170122
#start.index <- 270 # for 21:13 to 23:55 Z
#end.index   <- 597
# NOTE 234611.nc MIXPHA file is missing.  CRASHES ON INDEX NUM 576

# for 20170131
#start.index <- 90  # for 20:19 to ~22:10 Z
#end.index   <- 317

# for 20170204
#start.index <- 79  # for 22:08 to 22:30 Z
#end.index   <- 123
#start.index <- 182  # for 22:59 to 23:18 Z
#end.index   <- 219

# for 20170207
#start.index <- 83  # for 20:12 to 21:12 Z THIS IS ADAM'S TIME PERIOD
#end.index   <- 204
#start.index <- 205  # for 21:12 to 23:00 Z **** SHOULD INCLUDE THIS PERIOD? ****
#end.index   <- 420

# for 20170220
#start.index <- 58  # for 14:50 to 15:04 Z
#end.index   <- 86

# for 20170221
#start.index <- 49 # for 14:43 to 16:27 Z  
#end.index   <- 256

# for 20170309
#start.index <- 35  # for 14:03 to 14:15Z
#end.index   <- 55  # for 14:03 to 14:15Z
#start.index <- 143 # for 14:57 to 17Z
#end.index   <- 460 # 20:21 Z onward
#start.index <- 461 # for 14:57 to 17Z
#end.index   <- 886 # 20:21 Z onward

#----------------------------------------------
# define data paths
#----------------------------------------------
# define RadIA SRPC nc interest directories
radia.SRPC.nc.FZDZ.dir             <- file.path("/d1/serke/projects/case_studies/SNOWIE/data/RadIA_data/SRPC/nc/FRZDRZ/")
radia.SRPC.nc.MIXPHA.dir           <- file.path("/d1/serke/projects/case_studies/SNOWIE/data/RadIA_data/SRPC/nc/MIXPHA/")
radia.SRPC.nc.SLW.dir              <- file.path("/d1/serke/projects/case_studies/SNOWIE/data/RadIA_data/SRPC/nc/SLW/")
radia.SRPC.nc.PLATES.dir           <- file.path("/d1/serke/projects/case_studies/SNOWIE/data/RadIA_data/SRPC/nc/PLATES/")

# define 10/30s obs data file path
# NOTE: THESE DATA HAVE BEEN MOVED TO /media/serke/pens_d2/d1/serke/projects/case_studies/SNOWIE/data/AC_data/UWKA_data/30s_average
filename.30s.obs                   <- file.path("/d1/serke/projects/case_studies/SNOWIE/data/AC_data/UWKA_data/30s_average/v3.2/snowieDmax30sObservations.v3.txt")

# define Korolev ice/liq binned data file path
koro.2DS.sizedist.dir              <- file.path("/d1/serke/projects/case_studies/SNOWIE/data/AC_data/UWKA_data/2DS_part_probe/reprocessed_korolev/SNOWIE_2DS_SizeDistribution_v190501/")
koro.2DS.sizedist.ice.file         <- paste("2DS_SizeDistribution_Dmax_Ice_", radia.date, ".csv", sep="")
koro.2DS.sizedist.liq.file         <- paste("2DS_SizeDistribution_Dmax_Liq_", radia.date, ".csv", sep="")

# define NEXRAD site location csv file path
nexrad.site.dataframetxt.dir       <- file.path("/d1/serke/projects/case_studies/SNOWIE/data/RadIA_data/nexrad_site_data/")

# output data
output.file.dir                    <- file.path("/d1/serke/projects/case_studies/SNOWIE/data/RadIA_matchto_AC_data/v3.2/")

#----------------------------------------------
# define a list of nc files in dir
#----------------------------------------------
radia.SRPC.nc.FRZDZ.filelist.str   <- list.files(paste(radia.SRPC.nc.FZDZ.dir, radia.date, sep = ""),  pattern = "[.nc]",       full.names = FALSE, ignore.case = TRUE)
radia.SRPC.nc.FRZDZ.hhmmsslist.str <- substring(radia.SRPC.nc.FRZDZ.filelist.str, 1, 6)
radia.SRPC.nc.FRZDZ.hhmmsslist.num <- unlist(lapply(radia.SRPC.nc.FRZDZ.hhmmsslist.str, as.numeric))

# print results to screen for debugging
radia.SRPC.nc.FRZDZ.filelist.str 
#radia.SRPC.nc.FRZDZ.hhmmsslist.str
#radia.SRPC.nc.FRZDZ.hhmmsslist.num

#----------------------------------------------
# load the NEXRAD site location text file
#----------------------------------------------
NEXRAD.site.df                    <- read.csv(paste(nexrad.site.dataframetxt.dir, "nexrad_site.csv", sep = ""), header = FALSE, sep = ",", dec = ".", stringsAsFactors=FALSE)
colnames(NEXRAD.site.df)          <- c("NCDCID", "ICAO", "WBAN", "radname", "COUNTRY", "STATE", "COUNTY", "lat", "lon", "elev", "GMTdiff", "STN_TYPE")
head(NEXRAD.site.df)

# load info specific to KCBX
ind.KCBX                          <- which(NEXRAD.site.df$ICAO == " KCBX")
KCBX.elev.ft                      <- NEXRAD.site.df$elev[ind.KCBX]
KCBX.lat.deg                      <- NEXRAD.site.df$lat[ind.KCBX]
KCBX.lon.deg                      <- NEXRAD.site.df$lon[ind.KCBX]

#----------------------------------------------
# load 10 or 30 second averaged observational data set into data frame
#----------------------------------------------
# for 30s obs data, version 2 (w/ koro values)
rawData                            <- read.csv(filename.30s.obs, header = TRUE, sep = ",", dec = ".") 
snowie.Dmax.30s.obs.colnames       <- c(           "lat.deg",   "lon.deg",         "alt.m",       "time",         "d99.liq.gt50.um",         "d99.liq.gt10.um",         "dmax.liq.um",          "d99.ice.gt50.um",         "d99.ice.gt10.um",         "dmax.ice.um",                  "2DS.num.liq.gt50",                 "2DS.num.liq.gt10",                 "2DS.num.ice.gt50",                 "2DS.num.ice.gt10",        "num.liq.frac.gt50",       "num.liq.frac.gt10",                "2DS.mass.liq.gt50",                "2DS.mass.liq.gt10",                "2DS.mass.ice.gt50",                "2DS.mass.ice.gt10",       "mass.liq.frac.gt50",       "mass.liq.frac.gt10",                "cdp.num",                  "head.deg",      "revT.degC",      "rosT.degC",           "lwc100.gm3",           "lwcrose.gm3",           "lwcpvm.gm3",           "lwccdp.gm3",           "lwcnev.gm3",           "twcnev.gm3",          "meanvv.ms",           "rhlicor",           "rh",      "dewlicor.degC",      "dew.degC",             "pres.hpa")
snowie.Dmax.30s.obs.df             <- data.frame(rawData$lat, rawData$lon, rawData$alt..m., rawData$time, rawData$d99.liq.gt50..um., rawData$d99.liq.gt10..um.,  rawData$dmax.liq..um., rawData$d99.ice.gt50..um., rawData$d99.ice.gt10..um., rawData$dmax.ice..um.,  rawData$X2DS.num.liq.gt50....cm.3., rawData$X2DS.num.liq.gt10....cm.3., rawData$X2DS.num.ice.gt50....cm.3., rawData$X2DS.num.ice.gt10....cm.3.,  rawData$num.liq.frac.gt50, rawData$num.liq.frac.gt10, rawData$X2DS.mass.liq.gt50..g.m.3., rawData$X2DS.mass.liq.gt10..g.m.3., rawData$X2DS.mass.ice.gt50..g.m.3., rawData$X2DS.mass.ice.gt10..g.m.3., rawData$mass.liq.frac.gt50, rawData$mass.liq.frac.gt10,  rawData$cdpNum....cm.3.,  rawData$heading..deg.true., rawData$revT..C., rawData$rosT..C., rawData$lwc100..g.m.3., rawData$lwcRose..g.m.3., rawData$lwcPvm..g.m.3., rawData$lwcCdp..g.m.3., rawData$lwcNev..g.m.3., rawData$twcNev..g.m.3., rawData$meanvv..m.s., rawData$rhLicor...., rawData$rh...., rawData$dewLicor..C., rawData$dew..C., rawData$pressure..hPa.)
colnames(snowie.Dmax.30s.obs.df)   <- snowie.Dmax.30s.obs.colnames

# split time field into date and time components
yyyymmdd.str                       <- substring(as.character(rawData$time), 1, 11)
yymmdd.str                         <- format(strptime(yyyymmdd.str, format = "%d-%b-%Y"), "%Y-%m-%d")
hhmmss.str                         <- substring(as.character(rawData$time), 13, 20)

# make date and time components with a space between, as needed by 'asPOSIXct' function
yymmddhhmmss.wspace.str          <- paste(yymmdd.str, " ", hhmmss.str, sep="")

# add DateTime to data frame
snowie.Dmax.30s.obs.df$DateTime    <- yymmddhhmmss.wspace.str
snowie.Dmax.30s.obs.df$month       <- as.numeric(substring(yymmdd.str, 6, 7))
snowie.Dmax.30s.obs.df$day         <- as.numeric(substring(yymmdd.str, 9, 10))
#snowie.Dmax.10s.obs.df$DateTime    <- yymmddhhmmss.wspace.str

# view the data frame contents
head(snowie.Dmax.30s.obs.df)
#head(snowie.Dmax.10s.obs.df)

# search data frame for the given mmdd date
ind.mmdd                           <- which(snowie.Dmax.30s.obs.df$month == as.numeric(radia.mm) & snowie.Dmax.30s.obs.df$day == as.numeric(radia.dd))
num.ind.mmdd                       <- length(ind.mmdd)
print(paste(num.ind.mmdd, " data points on [mmdd] ", radia.date, sep=""))

#----------------------------------------------
# load 2DS_SizeDist ice and liq data into data frame
#----------------------------------------------
#koro.2DS.sizedist.ice.df          <- read.csv(paste(koro.2DS.sizedist.dir, koro.2DS.sizedist.ice.file, sep = ""), header = TRUE, sep = ",", dec = ".", stringsAsFactors=FALSE)
#koro.2DS.sizedist.liq.df          <- read.csv(paste(koro.2DS.sizedist.dir, koro.2DS.sizedist.liq.file, sep = ""), header = TRUE, sep = ",", dec = ".", stringsAsFactors=FALSE)

## convert seconds since start of day to POSIX
#koro.2DS.sizedist.ice.df$POSIX    <- as.POSIXct(trunc(Sys.time(), units="days") + koro.2DS.sizedist.liq.df$time.s.[1:8694])
#koro.2DS.sizedist.liq.df$POSIX    <- as.POSIXct(trunc(Sys.time(), units="days") + koro.2DS.sizedist.liq.df$time.s.)

#drops                             <- c("time.s.", "POSIX")
#koro.2DS.sizedist.liq.onlybins.df <- koro.2DS.sizedist.liq.df[ , !(names(koro.2DS.sizedist.liq.df) %in% drops)]
#koro.2DS.sizedist.ice.onlybins.df <- koro.2DS.sizedist.ice.df[ , !(names(koro.2DS.sizedist.ice.df) %in% drops)]

##head(koro.2DS.sizedist.ice.df)

#bin.center.micron                 <- seq(10, 2560, by=10)

lastncfile                        <- "200000.nc"

#----------------------------------------------
# loop through all of the 10s obs data lines of ind.mmdd for the given date to find the hh:mm:ss values and match them to hh mm ss in the filelist of nc files
#----------------------------------------------
#for (k in 1:num.ind.mmdd) {
#for (k in start.index:num.ind.mmdd) { 
for (k in start.index:end.index) {

  # index each line from desired date in 10s.obs df to exact koro profile
  #k  <-  100
  radia.hhmmss.str                   <- substring(snowie.Dmax.30s.obs.df$DateTime[ind.mmdd[k]], 12, 19)
  radia.hh.num                       <- as.numeric(substring(radia.hhmmss.str, 1, 2))
  radia.mm.num                       <- as.numeric(substring(radia.hhmmss.str, 4, 5))
  radia.ss.num                       <- as.numeric(substring(radia.hhmmss.str, 7, 8))
  radia.hhmmss.num                   <- as.numeric(paste(substring(radia.hhmmss.str, 1, 2), substring(radia.hhmmss.str, 4, 5), substring(radia.hhmmss.str, 7, 8), sep=""))
  radia.hhmmss.num
  #ind.time.match.koro                <- which(substring(koro.2DS.sizedist.liq.df$POSIX, 12, 19) == radia.hhmmss.str )
  
  #
  # which nc file has the minimum hhmmss time difference
  ind.closestncfile.to.30s.obs       <- which.min(abs(radia.SRPC.nc.FRZDZ.hhmmsslist.num-radia.hhmmss.num))
  ind.closestncfile.to.30s.obs

  # here is the string name of the closest nc file in time
  closestncfile.to.30s.obs           <- radia.SRPC.nc.FRZDZ.filelist.str[ind.closestncfile.to.30s.obs]

  # here is the time difference, in hhmmss
  timediff.closestncfile.to.30s.obs  <- abs(radia.SRPC.nc.FRZDZ.hhmmsslist.num[ind.closestncfile.to.30s.obs] - radia.hhmmss.num)
  print("---------------------")
  print(paste(" Time difference between ", as.character(radia.date), "'s ", as.character(k), "th 30-sec obs at time ", radia.hhmmss.str,  " and the closest RadIA nc file is ", as.character(timediff.closestncfile.to.30s.obs), " [hhmmss]", sep=""))
  
  # skip loading of RadIA INT nc files and all cartesian lat/lon calcs because these fields all exist from the last matchup time
  if (closestncfile.to.30s.obs == lastncfile) {
    print("closest == last, skip loading and cartesian, lat, lon calculations")
    
  } else {
    
    #----------------------------------------------
    # load RadIA interest files 
    #----------------------------------------------
    # FZDZ
    nc.radia.frzdrz.filename           <- paste(radia.SRPC.nc.FZDZ.dir, radia.date, "/", closestncfile.to.30s.obs, sep="")
    print(paste("loading ", nc.radia.frzdrz.filename, sep=""))
    if (substr(nc.radia.frzdrz.filename, nchar(nc.radia.frzdrz.filename)-1, nchar(nc.radia.frzdrz.filename)) == "gz") {
      print("File is gzipped. Unzipping...")
      untar(nc.radia.frzdrz.filename)
      nc.radia.frzdrz.filename <- substr(nc.radia.frzdrz.filename, 1, nchar(nc.radia.frzdrz.filename)-3)
    }
    nc.radia.frzdrz                    <- nc_open(nc.radia.frzdrz.filename, write = FALSE, verbose = FALSE)
    print(paste("The file has", nc.radia.frzdrz$nvars, "variables"))
    radia.frzdrz.var.num               <- seq(1, nc.radia.frzdrz$nvars, by=1)
    for (i in 1:length(radia.frzdrz.var.num)) {
      radia.frzdrz.nam <- paste("v", radia.frzdrz.var.num[i], sep = "")
      assign(radia.frzdrz.nam, nc.radia.frzdrz$var[[radia.frzdrz.var.num[i]]])
    }
    if (nc.radia.frzdrz$nvars == 9) {
      MeanDBZ                            <- ncvar_get( nc.radia.frzdrz, v1 )
      FRZDRZ                             <- ncvar_get( nc.radia.frzdrz, v9 )
    } else if (nc.radia.frzdrz$nvars == 1) {
      FRZDRZ                             <- ncvar_get( nc.radia.frzdrz, v1 )
    }
    #print(paste("V1 has name", v1$name))
    nc_close(nc.radia.frzdrz)
  
    # temp not working: 20170309 1426z file didnt have SLW field in the var list...... diagnose and fix
    # SLW
    nc.radia.slw.filename           <- paste(radia.SRPC.nc.SLW.dir, radia.date, "/", closestncfile.to.30s.obs, sep="")
    print(paste("loading ", nc.radia.slw.filename, sep=""))
    if (substr(nc.radia.slw.filename, nchar(nc.radia.slw.filename)-1, nchar(nc.radia.slw.filename)) == "gz") {
      print("File is gzipped. Unzipping...")
      untar(nc.radia.slw.filename)
      nc.radia.slw.filename <- substr(nc.radia.slw.filename, 1, nchar(nc.radia.slw.filename)-3)
    }
    nc.radia.slw                    <- nc_open(nc.radia.slw.filename, write = FALSE, verbose = FALSE)
    print(paste("The file has", nc.radia.slw$nvars, "variables"))
    radia.slw.var.num               <- seq(1, nc.radia.slw$nvars, by=1)
    for (i in 1:length(radia.slw.var.num)) {
      radia.slw.nam <- paste("v", radia.slw.var.num[i], sep = "")
      assign(radia.slw.nam, nc.radia.slw$var[[radia.slw.var.num[i]]])
    }
    if (nc.radia.slw$nvars == 9) {
      MeanZDR                         <- ncvar_get( nc.radia.slw, v1 )
      MeanKDP                         <- ncvar_get( nc.radia.slw, v4 )
      SLW                             <- ncvar_get( nc.radia.slw, v7 )
    } else if (nc.radia.slw$nvars == 1) {
      SLW                             <- ncvar_get( nc.radia.slw, v1 )
    }
    print(paste("V9 has name", v9$name))
    nc_close(nc.radia.slw)
  
    # MIXPHA
    nc.radia.mixpha.filename           <- paste(radia.SRPC.nc.MIXPHA.dir, radia.date, "/", closestncfile.to.30s.obs, sep="")
    print(paste("loading ", nc.radia.mixpha.filename, sep=""))
    if (substr(nc.radia.mixpha.filename, nchar(nc.radia.mixpha.filename)-1, nchar(nc.radia.mixpha.filename)) == "gz") {
      print("File is gzipped. Unzipping...")
      untar(nc.radia.mixpha.filename)
      nc.radia.mixpha.filename <- substr(nc.radia.mixpha.filename, 1, nchar(nc.radia.mixpha.filename)-3)
    }
    nc.radia.mixpha                    <- nc_open(nc.radia.mixpha.filename, write = FALSE, verbose = FALSE)
    print(paste("The file has", nc.radia.mixpha$nvars, "variables"))
    radia.mixpha.var.num               <- seq(1, nc.radia.mixpha$nvars, by=1)
    for (i in 1:length(radia.mixpha.var.num)) {
      radia.mixpha.nam <- paste("v", radia.mixpha.var.num[i], sep = "")
      assign(radia.mixpha.nam, nc.radia.mixpha$var[[radia.mixpha.var.num[i]]])
    }
    if (nc.radia.mixpha$nvars == 1) {
      MIXPHA                         <- ncvar_get( nc.radia.mixpha, v1 )
      print(paste("V1 has name", v1$name))
    } else if (nc.radia.mixpha$nvars == 4) {
      MIXPHA                         <- ncvar_get( nc.radia.mixpha, v4 )
      print(paste("V4 has name", v4$name))
    }
    nc_close(nc.radia.mixpha)
    
    #-------------------------------------------------
    # define polar coordinate fields: tilt angles (theta), azim and range data
    #-------------------------------------------------
    # define tilts (theta) in volume VCP
  
    theta.deg                          <- nc.radia.frzdrz$dim[2]$height$vals
    theta.rad                          <- NISTdegTOradian(theta.deg)
  
    # define ranges
    radia.range.atom.km                <- seq(2.125, 459.875, 0.25)            # length is 1832
    radia.range.matr.km                <- rep.col(radia.range.atom.km, 720)
  
    # define azimuths
    radia.azim.atom.deg                <- seq(0.0, 359.5, 0.5)              # length is 720
    radia.azim.matr.deg                <- rep.row(radia.azim.atom.deg, 1832)

    #-------------------------------------------------
    # manipulate polar coordinate fields
    #-------------------------------------------------
    # convert range to altitude for each given tilt angle (theta), accounting for curvature of the earth
    num.theta.deg                      <- length(nc.radia.frzdrz$dim[2]$height$vals)
    alt.beam.km                        <- replicate(num.theta.deg, radia.range.matr.km)   # make a matrix for every tilt
    for (j in 1:length(theta.rad)) {
      alt.beam.km[,,j]                 <- (radia.range.matr.km * sin(theta.rad[j])) + ((radia.range.matr.km^2)/ radius.earth.km)
    }

    # convert polar coordinates to cartesian coordinates
    radia.coords.cart                  <- pol2cart(as.vector(radia.range.matr.km), as.vector(NISTdegTOradian(radia.azim.matr.deg))) #pol2cart from 'useful'
    x                                  <- radia.coords.cart$x
    y                                  <- radia.coords.cart$y
    radia.coords                       <- data.frame(x, y)

    # compute lat/lon for every x/y pair in radia.coords.xy and add array to data frame
    ninetydeg.rad                      <- 90*pi/180
    lat.0.rad                          <- NISTdegTOradian(KCBX.lat.deg)
    lon.0.rad                          <- NISTdegTOradian(abs(KCBX.lon.deg))
    lat.new.deg                        <- rep(0, length(radia.coords$x))
    lon.new.deg                        <- rep(0, length(radia.coords$x))

    print("Into kk loop")
    for (kk in 1:length(radia.coords$x)) {
      anglesxy        <- atan(radia.coords$y[kk] / radia.coords$x[kk])
      tot.dist.km     <- radia.coords$x[kk]      / cos(anglesxy)         # in km from radar
      bearing.rad     <- anglesxy + ninetydeg.rad
      lat.new.rad     <- asin(sin(lat.0.rad) * cos(tot.dist.km/radius.earth.km) + cos(lat.0.rad) * sin(tot.dist.km/radius.earth.km) * cos(bearing.rad))
      lon.new.rad     <- lon.0.rad + atan2(sin(bearing.rad) * sin(tot.dist.km/radius.earth.km) * cos(lat.0.rad), cos(tot.dist.km/radius.earth.km) * sin(lat.new.rad) )
      y[kk]           <- NISTradianTOdeg(lat.new.rad)
      x[kk]           <- NISTradianTOdeg(lon.new.rad)
      #print(lat.new.deg[kk])
      #print(lon.new.deg[kk])
    }
  
  } # end of if (closestncfile.to.30s.obs == lastncfile) 

  lastncfile                     <- closestncfile.to.30s.obs
  
  # save computed lat/lon into radia.coords df
  radia.coords$lat                   <- as.data.frame(y)
  radia.coords$lon                   <- as.data.frame(-x)
  lat                                <- y
  lon                                <- x
  radia.coords$location.id           <- seq(1, length(x), 1)  # don't think this is needed anymore
  #colnames(radia.coords)             <- c(x, y, lat, lon)
  head(radia.coords)

  # calculate distance from each radar point to the radar
  #dist.geo.radia.pt.to.radar     <- distGeo(matrix(c(as.numeric(unlist(radia.coords$lon)), as.numeric(unlist(radia.coords$lat))), ncol=2), matrix(c(KCBX.lon.deg, KCBX.lat.deg), ncol=2)) # in meters
  # find the minimum distance
  #ind.dist.geo.radia.pt.to.radar <- which.min(dist.geo.radia.pt.to.radar)
  
  #------------------------------------------------------------
  # indexing to match radia domain to this kth aircraft values
  #------------------------------------------------------------
  # calculate distance from each radar point to the current aircraft location
  print("Into distGeo")
  dist.geo.radia.pt.to.ac        <- distGeo(matrix(c(as.numeric(unlist(radia.coords$lon)), as.numeric(unlist(radia.coords$lat))), ncol=2), matrix(c(snowie.Dmax.30s.obs.df$lon[ind.mmdd[k]], snowie.Dmax.30s.obs.df$lat[ind.mmdd[k]]), ncol=2)) # in meters
 
  # find the index of the RadIA points with minimum distance from the kth AC point 
  ind.dist.geo.radia.pt.to.ac    <- which.min(dist.geo.radia.pt.to.ac)

  # WE HAVE FOUND AND LOADED THE NC FILE THAT IS CLOSEST IN TIME TO THE AC POINT.  
  # WE HAVE CONVERTED THE RADAR POLAR VOLUME TO CARTESIAN COORDINATES. 
  # WE HAVE FOUND THE HEIGHT OF EVERY RADAR PIXEL ABOVE A CURVED EARTH SURFACE.  
  # WE HAVE FOUND THE SMALLEST DISTANCE BETWEEN AIRCRAFT POINT AND RADAR PIXEL.
  
  # NEXT STEP IS TO FIND THE CLOSEST RADAR COLUMN TO THE AIRCRAFT LOCATION, THEN FIND THE CLOSEST HEIGHT TO THE RADAR BEAM HEIGHT AND EXTRACT THE BEST REPRESENTATION OF THAT MATCHUP
  as.numeric(unlist(radia.coords$lon))[ind.dist.geo.radia.pt.to.ac] # RadIA coords closest to AC
  as.numeric(unlist(radia.coords$lat))[ind.dist.geo.radia.pt.to.ac]
  as.numeric(unlist(radia.coords$x))[ind.dist.geo.radia.pt.to.ac]
  as.numeric(unlist(radia.coords$y))[ind.dist.geo.radia.pt.to.ac]
  
  print("Into distm radar to AC")
  dist.radar.to.ac.m             <- distm(matrix(c(KCBX.lon.deg, KCBX.lat.deg), ncol=2), matrix(c(snowie.Dmax.30s.obs.df$lon[ind.mmdd[k]], snowie.Dmax.30s.obs.df$lat[ind.mmdd[k]]), ncol=2), fun = distHaversine)
  rangegates.radar.to.ac         <- ((dist.radar.to.ac.m/1000) - 0.5) / 0.25  # convert m to km, subtract near field idst (0.5) and divide by range resolution (0.25)
  azim.radar.to.ac.deg           <- bearing(c(KCBX.lon.deg, KCBX.lat.deg), c(snowie.Dmax.30s.obs.df$lon[ind.mmdd[k]], snowie.Dmax.30s.obs.df$lat[ind.mmdd[k]]))
  if (azim.radar.to.ac.deg < 0) {
    azim.radar.to.ac.deg  <- 360 + azim.radar.to.ac.deg
  }
  print(paste(" The AC is ", as.character(round(rangegates.radar.to.ac, digits=2)), " range gates away at ", as.character(round(azim.radar.to.ac.deg, digits=2)), " degrees azim from radar", sep=""))
  ind.closest.azim.deg           <- which.min(abs(radia.azim.atom.deg - azim.radar.to.ac.deg))
  if (ind.closest.azim.deg < 6) {
    ind.closest.azim.deg <- 6
  } else if (ind.closest.azim.deg > 714) {
    ind.closest.azim.deg <- 715
  }
  ind.closest.range.num          <- round(rangegates.radar.to.ac, digits=0)
  
  #------------------------------------------------------
  # INTERESTS: MAX_IN_COLUMN VALUES CLOSEST TO AC POSITION
  #------------------------------------------------------
  # FZDZ 
  FZDZ.int.match.column          <- FRZDRZ[ind.closest.range.num, ind.closest.azim.deg, ]
  if (all(is.na(FZDZ.int.match.column), TRUE) == FALSE) {
    FZDZ.int.match.column.max    <- round(max(FZDZ.int.match.column, na.rm=TRUE), digits=2)
  } else {
    FZDZ.int.match.column.max    <- NA
  }
  #print(paste("FZDZ-max-in-column closest to AC is: ", as.character(FZDZ.int.match.column.max), sep=""))
  
  # temp disabled
  # SLW 
  SLW.int.match.column           <- SLW[ind.closest.range.num, ind.closest.azim.deg, ]
  if (all(is.na(SLW.int.match.column), TRUE) == FALSE) {
    SLW.int.match.column.max    <- round(max(SLW.int.match.column, na.rm=TRUE), digits=2)
  } else {
    SLW.int.match.column.max    <- NA
  }
  #print(paste("SLW-max-in-column closest to AC is: ", as.character(SLW.int.match.column.max), sep=""))
  
  # MIXPHA 
  MIXPHA.int.match.column        <- MIXPHA[ind.closest.range.num, ind.closest.azim.deg, ]
  if (all(is.na(MIXPHA.int.match.column), TRUE) == FALSE) {
    MIXPHA.int.match.column.max    <- round(max(MIXPHA.int.match.column, na.rm=TRUE), digits=2)
  } else {
    MIXPHA.int.match.column.max    <- NA
  }
  #print(paste("MIXPHA-max-in-column closest to AC is: ", as.character(MIXPHA.int.match.column.max), sep=""))
  #-----------------------------------------------------
  
  #-----------------------------------------------------
  # INTERESTS: 11 x 11 COLUMN VALUES CENTERED ON CLOSEST TO AC COLUMN
  #-----------------------------------------------------
  # FZDZ
  FZDZ.int.match.matrix          <- FRZDRZ[(ind.closest.range.num-5):(ind.closest.range.num+5), (ind.closest.azim.deg-5):(ind.closest.azim.deg+5), ]
  FZDZ.int.match.matrix.max      <- round(max(FZDZ.int.match.matrix, na.rm=TRUE),  digits=2)
  FZDZ.int.match.matrix.mean     <- round(mean(FZDZ.int.match.matrix, na.rm=TRUE), digits=2)
  if (is.infinite(FZDZ.int.match.matrix.max)) {
    FZDZ.int.match.matrix.max <- NaN
  }
  if (is.infinite(FZDZ.int.match.matrix.mean)) {
    FZDZ.int.match.matrix.mean <- NaN
  }
  #print(paste("FZDZ-max11x11 closest to AC is: ", as.character(FZDZ.int.match.matrix.max), sep=""))
  #print(paste("FZDZ-mean11x11 closest to AC is: ", as.character(FZDZ.int.match.matrix.mean), sep=""))
  
  # temp disabled
  # SLW
  SLW.int.match.matrix          <- SLW[(ind.closest.range.num-5):(ind.closest.range.num+5), (ind.closest.azim.deg-5):(ind.closest.azim.deg+5), ]
  SLW.int.match.matrix.max      <- round(max(SLW.int.match.matrix, na.rm=TRUE),  digits=2)
  SLW.int.match.matrix.mean     <- round(mean(SLW.int.match.matrix, na.rm=TRUE), digits=2)
  if (is.infinite(SLW.int.match.matrix.max)) {
    SLW.int.match.matrix.max <- NaN
  }
  if (is.infinite(SLW.int.match.matrix.mean)) {
    SLW.int.match.matrix.mean <- NaN
  }
  #print(paste("SLW-max11x11 closest to AC is: ", as.character(SLW.int.match.matrix.max), sep=""))
  #print(paste("SLW-mean11x11 closest to AC is: ", as.character(SLW.int.match.matrix.mean), sep=""))
  
  # MIXPHA
  MIXPHA.int.match.matrix          <- MIXPHA[(ind.closest.range.num-5):(ind.closest.range.num+5), (ind.closest.azim.deg-5):(ind.closest.azim.deg+5), ]
  MIXPHA.int.match.matrix.max      <- round(max(MIXPHA.int.match.matrix, na.rm=TRUE),  digits=2)
  MIXPHA.int.match.matrix.mean     <- round(mean(MIXPHA.int.match.matrix, na.rm=TRUE), digits=2)
  if (is.infinite(MIXPHA.int.match.matrix.max)) {
    MIXPHA.int.match.matrix.max <- NaN
  }
  if (is.infinite(MIXPHA.int.match.matrix.mean)) {
    MIXPHA.int.match.matrix.mean <- NaN
  }
  #print(paste("MIXPHA-max11x11 closest to AC is: ", as.character(MIXPHA.int.match.matrix.max), sep=""))
  #print(paste("MIXPHA-mean11x11 closest to AC is: ", as.character(MIXPHA.int.match.matrix.mean), sep=""))
  #-----------------------------------------------------
  
  AC.alt.m                       <- round(snowie.Dmax.30s.obs.df$alt.m[ind.mmdd[k]],              digits=0)
  AC.lat.deg                     <- round(snowie.Dmax.30s.obs.df$lat.deg[ind.mmdd[k]],            digits=4)
  AC.lon.deg                     <- round(snowie.Dmax.30s.obs.df$lon.deg[ind.mmdd[k]],            digits=4)
  T.rev.C                        <- round(snowie.Dmax.30s.obs.df$revT.degC[ind.mmdd[k]],          digits=1)
  T.ros.C                        <- round(snowie.Dmax.30s.obs.df$rosT.degC[ind.mmdd[k]],          digits=1)
  LWC.ros                        <- round(snowie.Dmax.30s.obs.df$lwcrose.gm3[ind.mmdd[k]],        digits=2)
  LWC.nev                        <- round(snowie.Dmax.30s.obs.df$lwcnev.gm3[ind.mmdd[k]],         digits=2)
  LWC.cdp                        <- round(snowie.Dmax.30s.obs.df$lwccdp.gm3[ind.mmdd[k]],         digits=2)
  pres.mb                        <- round(snowie.Dmax.30s.obs.df$pres.hpa[ind.mmdd[k]],           digits=0)
  Dmax.um                        <- round(snowie.Dmax.30s.obs.df$dmax.liq.um[ind.mmdd[k]],        digits=0)
  if (is.nan(Dmax.um)) {
    Dmax.um                      <- 0
  } 
  D99.um                         <- round(snowie.Dmax.30s.obs.df$d99.liq.gt10[ind.mmdd[k]],       digits=0)
  if (is.nan(D99.um)) {
    D99.um                       <- 0
  } 
  mass.liq.frac                  <- round(snowie.Dmax.30s.obs.df$mass.liq.frac.gt10[ind.mmdd[k]], digits=3)
  if (is.null(snowie.Dmax.30s.obs.df$mass.liq.frac.gt10[ind.mmdd[k]])) {
    mass.liq.frac                <- -1
  } else if (is.nan(snowie.Dmax.30s.obs.df$mass.liq.frac.gt10[ind.mmdd[k]])) {
    mass.liq.frac                <- -1
  }
  
  #------------------------
  # CREATE DF FOR INTERESTS
  #------------------------
  # FZDZ fields
  #FZDZ.int.match.df              <- data.frame(FZDZ.int.match.column.max, FZDZ.int.match.matrix.max, FZDZ.int.match.matrix.mean)
  # SLW fields
  #SLW.int.match.df               <- data.frame(SLW.int.match.column.max, SLW.int.match.matrix.max, SLW.int.match.matrix.mean)
  # MIXPHA fields
  #MIXPHA.int.match.df            <- data.frame(MIXPHA.int.match.column.max, MIXPHA.int.match.matrix.max, MIXPHA.int.match.matrix.mean)
  
  # all fields
  if (k == start.index) {
    RadIA.ints.match.first.df      <- data.frame(radia.date, radia.hh.num, radia.mm.num, radia.ss.num, timediff.closestncfile.to.30s.obs, round(rangegates.radar.to.ac, digits=0), round(azim.radar.to.ac.deg, digits=2), AC.lat.deg, AC.lon.deg, AC.alt.m, T.rev.C, T.ros.C, pres.mb, LWC.ros, LWC.nev, LWC.cdp, mass.liq.frac, Dmax.um, D99.um, FZDZ.int.match.column.max, FZDZ.int.match.matrix.max, FZDZ.int.match.matrix.mean, SLW.int.match.column.max, SLW.int.match.matrix.max, SLW.int.match.matrix.mean, MIXPHA.int.match.column.max, MIXPHA.int.match.matrix.max, MIXPHA.int.match.matrix.mean)
    RadIA.ints.match.latest.df     <- data.frame(radia.date, radia.hh.num, radia.mm.num, radia.ss.num, timediff.closestncfile.to.30s.obs, round(rangegates.radar.to.ac, digits=0), round(azim.radar.to.ac.deg, digits=2), AC.lat.deg, AC.lon.deg, AC.alt.m, T.rev.C, T.ros.C, pres.mb, LWC.ros, LWC.nev, LWC.cdp, mass.liq.frac, Dmax.um, D99.um, FZDZ.int.match.column.max, FZDZ.int.match.matrix.max, FZDZ.int.match.matrix.mean, SLW.int.match.column.max, SLW.int.match.matrix.max, SLW.int.match.matrix.mean, MIXPHA.int.match.column.max, MIXPHA.int.match.matrix.max, MIXPHA.int.match.matrix.mean)                                                                                                                                                                             
    RadIA.ints.match.total.df      <- rbind(RadIA.ints.match.first.df, RadIA.ints.match.latest.df)
    rm(RadIA.ints.match.first.df)
  } else if (k > start.index) {
    RadIA.ints.match.latest.df     <- data.frame(radia.date, radia.hh.num, radia.mm.num, radia.ss.num, timediff.closestncfile.to.30s.obs, round(rangegates.radar.to.ac, digits=0), round(azim.radar.to.ac.deg, digits=2), AC.lat.deg, AC.lon.deg, AC.alt.m, T.rev.C, T.ros.C, pres.mb, LWC.ros, LWC.nev, LWC.cdp, mass.liq.frac, Dmax.um, D99.um, FZDZ.int.match.column.max, FZDZ.int.match.matrix.max, FZDZ.int.match.matrix.mean, SLW.int.match.column.max, SLW.int.match.matrix.max, SLW.int.match.matrix.mean, MIXPHA.int.match.column.max, MIXPHA.int.match.matrix.max, MIXPHA.int.match.matrix.mean)
    RadIA.ints.match.total.df      <- rbind(RadIA.ints.match.total.df, RadIA.ints.match.latest.df)
  }
  #if (k == 90) {
  #  RadIA.ints.match.first.df      <- data.frame(radia.date, radia.hh.num, radia.mm.num, radia.ss.num, timediff.closestncfile.to.30s.obs, round(rangegates.radar.to.ac, digits=0), round(azim.radar.to.ac.deg, digits=2), AC.lat.deg, AC.lon.deg, AC.alt.m, T.rev.C, T.ros.C, pres.mb, LWC.ros, LWC.nev, LWC.cdp, Dmax.um, D99.um, FZDZ.int.match.column.max, FZDZ.int.match.matrix.max, FZDZ.int.match.matrix.mean, SLW.int.match.column.max, SLW.int.match.matrix.max, SLW.int.match.matrix.mean, MIXPHA.int.match.column.max, MIXPHA.int.match.matrix.max, MIXPHA.int.match.matrix.mean)
  #  RadIA.ints.match.latest.df     <- data.frame(radia.date, radia.hh.num, radia.mm.num, radia.ss.num, timediff.closestncfile.to.30s.obs, round(rangegates.radar.to.ac, digits=0), round(azim.radar.to.ac.deg, digits=2), AC.lat.deg, AC.lon.deg, AC.alt.m, T.rev.C, T.ros.C, pres.mb, LWC.ros, LWC.nev, LWC.cdp, Dmax.um, D99.um, FZDZ.int.match.column.max, FZDZ.int.match.matrix.max, FZDZ.int.match.matrix.mean, SLW.int.match.column.max, SLW.int.match.matrix.max, SLW.int.match.matrix.mean, MIXPHA.int.match.column.max, MIXPHA.int.match.matrix.max, MIXPHA.int.match.matrix.mean)                                                                                                                                                                             
  #  RadIA.ints.match.total.df      <- rbind(RadIA.ints.match.first.df, RadIA.ints.match.latest.df)
  #  rm(RadIA.ints.match.first.df)
  #} else if (k > 90) {
  #  RadIA.ints.match.latest.df     <- data.frame(radia.date, radia.hh.num, radia.mm.num, radia.ss.num, timediff.closestncfile.to.30s.obs, round(rangegates.radar.to.ac, digits=0), round(azim.radar.to.ac.deg, digits=2), AC.lat.deg, AC.lon.deg, AC.alt.m, T.rev.C, T.ros.C, pres.mb, LWC.ros, LWC.nev, LWC.cdp, Dmax.um, D99.um, FZDZ.int.match.column.max, FZDZ.int.match.matrix.max, FZDZ.int.match.matrix.mean, SLW.int.match.column.max, SLW.int.match.matrix.max, SLW.int.match.matrix.mean, MIXPHA.int.match.column.max, MIXPHA.int.match.matrix.max, MIXPHA.int.match.matrix.mean)
  #  RadIA.ints.match.total.df      <- rbind(RadIA.ints.match.total.df, RadIA.ints.match.latest.df)
  #}
  rm(RadIA.ints.match.latest.df)
  #------------------------
  
  # write out data
  write.csv(RadIA.ints.match.total.df, file = paste(output.file.dir, radia.date, "_RadIA_ints_match.wSLW.v3.csv", sep=""), row.names=FALSE, na="")
  #write.csv(RadIA.ints.match.df, file = paste(output.file.dir, radia.date, "_", radia.hhmmss.num, "_RadIA_ints_match.csv", sep=""), row.names=TRUE, na="")
  
  # make an 11-second average around to center 10-sec obs time
  #ind.start                          <- ind.time.match.koro-5
  #ind.stop                           <- ind.time.match.koro+5
  #koro.2DS.sizedist.liq.onlybins.10s <- koro.2DS.sizedist.liq.onlybins.df[ind.start:ind.stop, ]
  #koro.2DS.sizedist.ice.onlybins.10s <- koro.2DS.sizedist.ice.onlybins.df[ind.start:ind.stop, 1:256]
  #for (i in 1: dim(koro.2DS.sizedist.liq.onlybins.10s)[2]) {
  #  koro.2DS.sizedist.liq.onlybins.10s.ave[i] <- mean(koro.2DS.sizedist.liq.onlybins.10s[, i])
  #  koro.2DS.sizedist.ice.onlybins.10s.ave[i] <- mean(koro.2DS.sizedist.ice.onlybins.10s[, i])
  #}
  
  # calculate Dmax for ice and liq
  #ind.thresh.liq <- which(koro.2DS.sizedist.liq.onlybins.10s.ave >= 0.1)
  #ind.thresh.ice <- which(koro.2DS.sizedist.ice.onlybins.10s.ave >= 0.1)
  #if(length(ind.thresh.liq) > 0) {
  #  dmax.liq       <- max(ind.thresh.liq) * 10
  #} else {
  #  dmax.liq       <- NA
  #}
  #if(length(ind.thresh.ice) > 0) {
  #  dmax.ice       <- max(ind.thresh.ice) * 10
  #} else {
  #  dmax.ice       <- NA
  #}
  
  # plot 10-s distribs of liq and ice conc from koro files
  #plot.new()
  #plot(bin.center.micron[1:25], koro.2DS.sizedist.liq.onlybins.10s.ave[1:25], type="o", col="red", main=paste(snowie.Dmax.10s.obs.df$DateTime[ind.mmdd[k]], " 10-sec mean", sep=""), xlab="bin center [um]", ylab="# / L")
  #if(koro.2DS.sizedist.ice.onlybins.10s.ave[1] > 10000) {
  #  koro.2DS.sizedist.ice.onlybins.10s.ave[1] = 0
  #}
  #lines(bin.center.micron[1:25], koro.2DS.sizedist.ice.onlybins.10s.ave[1:25], type="o", col="blue")
  #legend(150, max(koro.2DS.sizedist.liq.onlybins.df[ind.time.match.koro,1:25]*0.65), c("Liquid", "Solid"), col=c(2,4), pch=c(1,1))
  #text(  170, max(koro.2DS.sizedist.liq.onlybins.df[ind.time.match.koro,1:25]*0.51), paste("2DS Dmax-liq = ", dmax.liq,                                                         " um",    sep=""))
  #text(  170, max(koro.2DS.sizedist.liq.onlybins.df[ind.time.match.koro,1:25]*0.48), paste("2DS Dmax-ice = ", dmax.ice,                                                         " um",    sep=""))
  #text(  170, max(koro.2DS.sizedist.liq.onlybins.df[ind.time.match.koro,1:25]*0.45), paste("d99QC= ",         round(snowie.Dmax.10s.obs.df$d99QC.um[ind.mmdd[k]],    digits=0), " um",    sep=""))
  #text(  170, max(koro.2DS.sizedist.liq.onlybins.df[ind.time.match.koro,1:25]*0.42), paste("Frac-liq= ",      round(snowie.Dmax.10s.obs.df$fracliq[ind.mmdd[k]],     digits=2),           sep=""))
  #text(  170, max(koro.2DS.sizedist.liq.onlybins.df[ind.time.match.koro,1:25]*0.39), paste("Nevz LWC= ",      round(snowie.Dmax.10s.obs.df$lwcnev.gm3[ind.mmdd[k]],  digits=2), " gm-3",  sep=""))
  #text(  170, max(koro.2DS.sizedist.liq.onlybins.df[ind.time.match.koro,1:25]*0.36), paste("Rose LWC= ",      round(snowie.Dmax.10s.obs.df$lwcrose.gm3[ind.mmdd[k]], digits=2), " gm-3",  sep=""))
  #text(  170, max(koro.2DS.sizedist.liq.onlybins.df[ind.time.match.koro,1:25]*0.33), paste("FZDZ-clo-MIC= ",  round(FZDZ.int.match.column.max,                       digits=2),           sep=""))
  #text(  170, max(koro.2DS.sizedist.liq.onlybins.df[ind.time.match.koro,1:25]*0.30), paste("FZDZ-11^2-MIA= ", round(FZDZ.int.match.matrix.max,                       digits=2),           sep=""))
  #text(  170, max(koro.2DS.sizedist.liq.onlybins.df[ind.time.match.koro,1:25]*0.27), paste("FZDZ-11^2-AIA= ", round(FZDZ.int.match.matrix.mean,                      digits=2),           sep=""))
  #text(  170, max(koro.2DS.sizedist.liq.onlybins.df[ind.time.match.koro,1:25]*0.24), paste("FZSD-clo-MIC= ",  round(SLW.int.match.column.max,                        digits=2),           sep=""))
  #text(  170, max(koro.2DS.sizedist.liq.onlybins.df[ind.time.match.koro,1:25]*0.21), paste("FZSD-11^2-MIA= ", round(SLW.int.match.matrix.max,                        digits=2),           sep=""))
  #text(  170, max(koro.2DS.sizedist.liq.onlybins.df[ind.time.match.koro,1:25]*0.18), paste("FZSD-11^2-AIA= ", round(SLW.int.match.matrix.mean,                       digits=2),           sep=""))
  #text(  170, max(koro.2DS.sizedist.liq.onlybins.df[ind.time.match.koro,1:25]*0.15), paste("MPHA-clo-MIC= ",  round(MIXPHA.int.match.column.max,                     digits=2),           sep=""))
  #text(  170, max(koro.2DS.sizedist.liq.onlybins.df[ind.time.match.koro,1:25]*0.12), paste("MPHA-11^2-MIA= ", round(MIXPHA.int.match.matrix.max,                     digits=2),           sep=""))
  #text(  170, max(koro.2DS.sizedist.liq.onlybins.df[ind.time.match.koro,1:25]*0.09), paste("MPHA-11^2-AIA= ", round(MIXPHA.int.match.matrix.mean,                    digits=2),           sep=""))
  #text(  170, max(koro.2DS.sizedist.liq.onlybins.df[ind.time.match.koro,1:25]*0.06), paste("T-rev= ",         round(snowie.Dmax.10s.obs.df$revT.degC[ind.mmdd[k]],   digits=1), " deg C", sep=""))
  #text(  170, max(koro.2DS.sizedist.liq.onlybins.df[ind.time.match.koro,1:25]*0.03), paste("T-rev= ",         round(snowie.Dmax.10s.obs.df$rosT.degC[ind.mmdd[k]],   digits=1), " deg C", sep=""))
  #grid()
  
  
}  # end of for(k in 1: num.ind.mmdd)





