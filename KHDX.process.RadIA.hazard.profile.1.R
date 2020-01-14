#----------------------------------------------
#
# Name:    KHDX.process.RadIA.hazard.profile.1.R
#
# Purpose: 1) load radia INT files from case study
#          2) process profile of hazards over set location 
#
# Created: 6.20.2019 dserke
#
# Note:    Start working mostly in m, convert to kft for plotting
#----------------------------------------------

# add libraries
library(lubridate)
library(R.matlab)
library(ncdf4)
library(NISTunits)
library(useful)
library(gridExtra)

require(dplyr)
require(geosphere)

# remove scientific notation in printing
options(scipen=999)

#----------------------------------------------
# constants and parameters
#----------------------------------------------
# define yyyymmdd character string to work from
radia.yyyy.2                         <- "2019"
radia.mm.2                           <- "02"
radia.dd.2                           <- "23"
radia.date.2                         <- paste(radia.yyyy.2, radia.mm.2, radia.dd.2, sep = "")

radius.earth.km                      <- 6378.1

# conversion values
mileTOkft                            <-    5.28

# define left and right values of x in rectangular hazard profile plots
xleft                                <- 0
xright                               <- 1

# lat/lon/alt for three locations on WSMR
# sta=stallion aaf 
sta.lat.deg                          <-   33.8172   # degrees N
sta.lon.deg                          <- -106.6453   # degrees W
sta.alt.m                            <- 1501.0000   # m
sta.alt.ft                           <- NISTmeterTOft(sta.alt.m)

# ski=skillet knob
ski.lat.deg                          <-   33.2304   # degrees N
ski.lon.deg                          <- -106.6411   # degrees W
ski.alt.m                            <- 1432.0000   # m
ski.alt.ft                           <- NISTmeterTOft(ski.alt.m)

#  wst=white sands town
wst.lat.deg                          <-   32.3824   # degrees N
wst.lon.deg                          <- -106.4907   # degrees W
wst.alt.m                            <- 1300.0000   # m
wst.alt.ft                           <- NISTmeterTOft(wst.alt.m)


#----------------------------------------------
# get user input regarding number of neighboring points to use 
#----------------------------------------------
# number of ranges around profile point
user.num.ranges  <- readline(prompt="Enter number of adjacent ranges: ")
# number of azimuths around profile point
user.num.azims   <- readline(prompt="Enter number of adjacent azimuths: ")
# stats method to represent profile
#  ex: mean, median, maximum
user.stat.method <- readline(prompt="Enter stat method [ex: 'mean', 'median', 'max']: ")

#----------------------------------------------
# define data paths
#----------------------------------------------
# define RadIA SRPC nc interest directories
campaign                           <- "ATEC_ranges_2019"
radia.SRPC.nc.FZDZ.dir             <- file.path(paste("/d1/serke/projects/case_studies/", campaign, "/data/RadIA_data/SRPC/nc/FRZDRZ/", radia.date.2, sep=""))
radia.SRPC.nc.MIXPHA.dir           <- file.path(paste("/d1/serke/projects/case_studies/", campaign, "/data/RadIA_data/SRPC/nc/MIXPHA/", radia.date.2, sep=""))
radia.SRPC.nc.SLW.dir              <- file.path(paste("/d1/serke/projects/case_studies/", campaign, "/data/RadIA_data/SRPC/nc/SLW/",    radia.date.2, sep=""))
#radia.SRPC.nc.PLATES.dir           <- file.path(paste("/d1/serke/projects/case_studies/", campaign, "/data/RadIA_data/SRPC/nc/PLATES/", radia.date.2, sep=""))

# define NEXRAD site location csv file path
nexrad.site.dataframetxt.dir       <- file.path("/d1/serke/projects/case_studies/SNOWIE/data/RadIA_data/nexrad_site_data/")

# output data
output.file.dir                    <- file.path(paste("/d1/serke/projects/case_studies/", campaign, "/data/RadIA_data/SRPC/output"))

#----------------------------------------------
# define a list of nc files in dir
#----------------------------------------------
radia.SRPC.nc.FRZDZ.filelist.str   <- list.files(radia.SRPC.nc.FZDZ.dir,  pattern = "[.nc]",       full.names = FALSE, ignore.case = TRUE)
radia.SRPC.nc.FRZDZ.hhmmsslist.str <- substring(radia.SRPC.nc.FRZDZ.filelist.str, 1, 6)
radia.SRPC.nc.FRZDZ.hhmmsslist.num <- unlist(lapply(radia.SRPC.nc.FRZDZ.hhmmsslist.str, as.numeric))

# print results to screen for debugging
radia.SRPC.nc.FRZDZ.filelist.str 
#radia.SRPC.nc.FRZDZ.hhmmsslist.str
#radia.SRPC.nc.FRZDZ.hhmmsslist.num

#----------------------------------------------
# load the NEXRAD site location text file
#----------------------------------------------
NEXRAD.site.df                     <- read.csv(paste(nexrad.site.dataframetxt.dir, "nexrad_site.csv", sep = ""), header = FALSE, sep = ",", dec = ".", stringsAsFactors=FALSE)
colnames(NEXRAD.site.df)           <- c("NCDCID", "ICAO", "WBAN", "radname", "COUNTRY", "STATE", "COUNTY", "lat", "lon", "elev", "GMTdiff", "STN_TYPE")
head(NEXRAD.site.df)

# load info specific to KCBX
ind.KHDX                           <- which(NEXRAD.site.df$ICAO == " KHDX")
KHDX.elev.ft                       <- NEXRAD.site.df$elev[ind.KHDX]
KHDX.elev.m                        <- NISTftTOmeter(KHDX.elev.ft)
KHDX.lat.deg                       <- NEXRAD.site.df$lat[ind.KHDX]
KHDX.lon.deg                       <- NEXRAD.site.df$lon[ind.KHDX]

lastncfile                         <- "200000.nc"

# initialize these matrices for later usage
profile.FZDZ.at.sta                <- matrix(, nrow=length(radia.SRPC.nc.FRZDZ.filelist.str), ncol=9)
profile.SLW.at.sta                 <- profile.FZDZ.at.sta
profile.MIXPHA.at.sta              <- profile.FZDZ.at.sta
profile.FZDZ.at.ski                <- profile.FZDZ.at.sta
profile.SLW.at.ski                 <- profile.FZDZ.at.sta
profile.MIXPHA.at.ski              <- profile.FZDZ.at.sta
profile.FZDZ.at.wst                <- profile.FZDZ.at.sta
profile.SLW.at.wst                 <- profile.FZDZ.at.sta
profile.MIXPHA.at.wst              <- profile.FZDZ.at.sta

#----------------------------------------------
# loop through each available nc file time
#----------------------------------------------
for (k in 1:length(radia.SRPC.nc.FRZDZ.filelist.str)) {

  print(paste("k= ", k, " out of ", length(radia.SRPC.nc.FRZDZ.filelist.str), sep=""))
  
  #----------------------------------------------
  # load RadIA interest files 
  #----------------------------------------------
  # FZDZ
  nc.radia.frzdrz.filename           <- paste(radia.SRPC.nc.FZDZ.dir, "/", radia.SRPC.nc.FRZDZ.filelist.str[k], sep="")
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
  
  # SLW
  nc.radia.slw.filename           <- paste(radia.SRPC.nc.SLW.dir, "/", radia.SRPC.nc.FRZDZ.filelist.str[k], sep="")
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
  print(paste("V1 has name", v1$name))
  nc_close(nc.radia.slw)
  
  # MIXPHA
  nc.radia.mixpha.filename           <- paste(radia.SRPC.nc.MIXPHA.dir, "/", radia.SRPC.nc.FRZDZ.filelist.str[k], sep="")
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
  # define tilts (theta), retrieve tilts from volume VCP
  theta.deg                          <- nc.radia.frzdrz$dim[2]$height$vals
  theta.rad                          <- NISTdegTOradian(theta.deg)

  # define ranges
  radia.range.atom.km                <- seq(2.125, 459.875, 0.25)            # length is 1832
  radia.range.matr.km                <- rep.col(radia.range.atom.km, 720)
  
  #-------------------------------------------------
  # manipulate polar coordinate fields
  #-------------------------------------------------
  # convert range to altitude for each given tilt angle (theta), accounting for curvature of the earth
  num.theta.deg                      <- length(nc.radia.frzdrz$dim[2]$height$vals)
  alt.beam.km                        <- replicate(num.theta.deg, radia.range.matr.km)   # make a matrix for every tilt dim=(1832, 720, 9)
  for (j in 1:length(theta.rad)) {
    alt.beam.km[,,j]                  <- (radia.range.matr.km  * sin(theta.rad[j])) + ((radia.range.matr.km^2)/ radius.earth.km)
  }
  
  lastncfile                     <- radia.SRPC.nc.FRZDZ.filelist.str[k]

  #-----------------------------------------------------------
  # work with profiles over key points on WSMR
  #-----------------------------------------------------------
  # calculate distances from key locations on WSMR to the Holoman radar (KHDX)
  #dist.geo.sta.to.radar          <- NISTmeterTOft(distGeo(matrix(c(sta.lon.deg, sta.lat.deg), ncol=2), matrix(c(KHDX.lon.deg, KHDX.lat.deg), ncol=2))) / 1000 # in kft
  #dist.geo.wst.to.radar          <- NISTmeterTOft(distGeo(matrix(c(wst.lon.deg, wst.lat.deg), ncol=2), matrix(c(KHDX.lon.deg, KHDX.lat.deg), ncol=2))) / 1000 # in kft
  #dist.geo.ski.to.radar          <- NISTmeterTOft(distGeo(matrix(c(ski.lon.deg, ski.lat.deg), ncol=2), matrix(c(KHDX.lon.deg, KHDX.lat.deg), ncol=2))) / 1000 # in kft
  dist.geo.sta.to.radar          <- distGeo(matrix(c(sta.lon.deg, sta.lat.deg), ncol=2), matrix(c(KHDX.lon.deg, KHDX.lat.deg), ncol=2)) # in m
  dist.geo.wst.to.radar          <- distGeo(matrix(c(wst.lon.deg, wst.lat.deg), ncol=2), matrix(c(KHDX.lon.deg, KHDX.lat.deg), ncol=2)) # in m
  dist.geo.ski.to.radar          <- distGeo(matrix(c(ski.lon.deg, ski.lat.deg), ncol=2), matrix(c(KHDX.lon.deg, KHDX.lat.deg), ncol=2)) # in m
  
  # calculate headings from key locations on WSMR to the Holoman radar (KHDX)
  azim.sta.to.radar              <- 360 + bearing(c(KHDX.lon.deg, KHDX.lat.deg), c(sta.lon.deg, sta.lat.deg)) # in degress from KHDX
  azim.ski.to.radar              <- 360 + bearing(c(KHDX.lon.deg, KHDX.lat.deg), c(ski.lon.deg, ski.lat.deg)) # in degress from KHDX
  azim.wst.to.radar              <- 360 + bearing(c(KHDX.lon.deg, KHDX.lat.deg), c(wst.lon.deg, wst.lat.deg)) # in degrees from KHDX
  
  # indexing to match radia domain to key locations on WSMR
  ind.range.sta.to.radar         <- round(dist.geo.sta.to.radar / 250, digits=0)
  ind.azim.sta.to.radar          <- round(azim.sta.to.radar, digits=0) * 2        # radar azim 720 values for 360 degrees, 0.5 degree resolution
  profile.SLW.at.sta[k,]         <- SLW[   ind.range.sta.to.radar, ind.azim.sta.to.radar, ]
  profile.FZDZ.at.sta[k,]        <- FRZDRZ[ind.range.sta.to.radar, ind.azim.sta.to.radar, ]
  profile.MIXPHA.at.sta[k,]      <- MIXPHA[ind.range.sta.to.radar, ind.azim.sta.to.radar, ]
  
  ind.range.ski.to.radar         <- round(dist.geo.ski.to.radar / 250, digits=0)
  ind.azim.ski.to.radar          <- round(azim.ski.to.radar, digits=0) * 2
  profile.SLW.at.ski[k,]         <- SLW[   ind.range.ski.to.radar, ind.azim.ski.to.radar, ]
  profile.FZDZ.at.ski[k,]        <- FRZDRZ[ind.range.ski.to.radar, ind.azim.ski.to.radar, ]
  profile.MIXPHA.at.ski[k,]      <- MIXPHA[ind.range.ski.to.radar, ind.azim.ski.to.radar, ]
     
  ind.range.wst.to.radar         <- round(dist.geo.wst.to.radar / 250, digits=0)
  ind.azim.wst.to.radar          <- round(azim.wst.to.radar, digits=0) * 2
  profile.SLW.at.wst[k,]         <- SLW[   ind.range.wst.to.radar, ind.azim.wst.to.radar, ]
  profile.FZDZ.at.wst[k,]        <- FRZDRZ[ind.range.wst.to.radar, ind.azim.wst.to.radar, ]
  profile.MIXPHA.at.wst[k,]      <- MIXPHA[ind.range.wst.to.radar, ind.azim.wst.to.radar, ]
  
  ## write out data
  #write.csv(RadIA.ints.match.total.df, file = paste(output.file.dir, radia.date, "_RadIA_ints_match.wSLW.v3.csv", sep=""), row.names=FALSE, na="")
  ##write.csv(RadIA.ints.match.df, file = paste(output.file.dir, radia.date, "_", radia.hhmmss.num, "_RadIA_ints_match.csv", sep=""), row.names=TRUE, na="")
  
}  # end of for(k in 1: length(radia.SRPC.nc.FZDZ.filelist.str))

# profile of radar beam center heights at key locations on WSMR
#prof.beam.height.sta      <- c(alt.beam.kft[ind.range.sta.to.radar,1,1], alt.beam.kft[ind.range.sta.to.radar,1,2], alt.beam.kft[ind.range.sta.to.radar,1,3], alt.beam.kft[ind.range.sta.to.radar,1,4], alt.beam.kft[ind.range.sta.to.radar,1,5], alt.beam.kft[ind.range.sta.to.radar,1,6], alt.beam.kft[ind.range.sta.to.radar,1,7], 11.3, 13.8)
#prof.beam.height.ski      <- c(alt.beam.kft[ind.range.ski.to.radar,1,1], alt.beam.kft[ind.range.ski.to.radar,1,2], alt.beam.kft[ind.range.ski.to.radar,1,3], alt.beam.kft[ind.range.ski.to.radar,1,4], alt.beam.kft[ind.range.ski.to.radar,1,5], alt.beam.kft[ind.range.ski.to.radar,1,6], alt.beam.kft[ind.range.ski.to.radar,1,7], 5.4, 6.8)
#prof.beam.height.wst      <- c(alt.beam.kft[ind.range.wst.to.radar,1,1], alt.beam.kft[ind.range.wst.to.radar,1,2], alt.beam.kft[ind.range.wst.to.radar,1,3], alt.beam.kft[ind.range.wst.to.radar,1,4], alt.beam.kft[ind.range.wst.to.radar,1,5], alt.beam.kft[ind.range.wst.to.radar,1,6], alt.beam.kft[ind.range.wst.to.radar,1,7], 9.7, 11.8)
prof.beam.height.sta      <- c(alt.beam.km[ind.range.sta.to.radar,1,1], alt.beam.km[ind.range.sta.to.radar,1,2], alt.beam.km[ind.range.sta.to.radar,1,3], alt.beam.km[ind.range.sta.to.radar,1,4], alt.beam.km[ind.range.sta.to.radar,1,5], alt.beam.km[ind.range.sta.to.radar,1,6], alt.beam.km[ind.range.sta.to.radar,1,7], 11.3, 13.8)
prof.beam.height.ski      <- c(alt.beam.km[ind.range.ski.to.radar,1,1], alt.beam.km[ind.range.ski.to.radar,1,2], alt.beam.km[ind.range.ski.to.radar,1,3], alt.beam.km[ind.range.ski.to.radar,1,4], alt.beam.km[ind.range.ski.to.radar,1,5], alt.beam.km[ind.range.ski.to.radar,1,6], alt.beam.km[ind.range.ski.to.radar,1,7], 5.4, 6.8)
prof.beam.height.wst      <- c(alt.beam.km[ind.range.wst.to.radar,1,1], alt.beam.km[ind.range.wst.to.radar,1,2], alt.beam.km[ind.range.wst.to.radar,1,3], alt.beam.km[ind.range.wst.to.radar,1,4], alt.beam.km[ind.range.wst.to.radar,1,5], alt.beam.km[ind.range.wst.to.radar,1,6], alt.beam.km[ind.range.wst.to.radar,1,7], 9.7, 11.8)

# profile of mean heights between radar beam center heights at key locations on WSMR
prof.betw.beam.height.sta <- c(prof.beam.height.sta[1]-0.3, (prof.beam.height.sta[1]+prof.beam.height.sta[2])/2, (prof.beam.height.sta[2]+prof.beam.height.sta[3])/2, (prof.beam.height.sta[3]+prof.beam.height.sta[4])/2, (prof.beam.height.sta[4]+prof.beam.height.sta[5])/2, (prof.beam.height.sta[5]+prof.beam.height.sta[6])/2, (prof.beam.height.sta[6]+prof.beam.height.sta[7])/2, (prof.beam.height.sta[7]+prof.beam.height.sta[8])/2, (prof.beam.height.sta[8]+prof.beam.height.sta[9])/2, prof.beam.height.sta[9]+0.3 )
prof.betw.beam.height.ski <- c(prof.beam.height.ski[1]-0.3, (prof.beam.height.ski[1]+prof.beam.height.ski[2])/2, (prof.beam.height.ski[2]+prof.beam.height.ski[3])/2, (prof.beam.height.ski[3]+prof.beam.height.ski[4])/2, (prof.beam.height.ski[4]+prof.beam.height.ski[5])/2, (prof.beam.height.ski[5]+prof.beam.height.ski[6])/2, (prof.beam.height.ski[6]+prof.beam.height.ski[7])/2, (prof.beam.height.ski[7]+prof.beam.height.ski[8])/2, (prof.beam.height.ski[8]+prof.beam.height.ski[9])/2, prof.beam.height.ski[9]+0.3)
prof.betw.beam.height.wst <- c(prof.beam.height.wst[1]-0.3, (prof.beam.height.wst[1]+prof.beam.height.wst[2])/2, (prof.beam.height.wst[2]+prof.beam.height.wst[3])/2, (prof.beam.height.wst[3]+prof.beam.height.wst[4])/2, (prof.beam.height.wst[4]+prof.beam.height.wst[5])/2, (prof.beam.height.wst[5]+prof.beam.height.wst[6])/2, (prof.beam.height.wst[6]+prof.beam.height.wst[7])/2, (prof.beam.height.wst[7]+prof.beam.height.wst[8])/2, (prof.beam.height.wst[8]+prof.beam.height.wst[9])/2, prof.beam.height.wst[9]+0.3)

# build out the data frame
sta.df                    <- data.frame(profile.FZDZ.at.sta, profile.SLW.at.sta, profile.MIXPHA.at.sta)
ski.df                    <- data.frame(profile.FZDZ.at.ski, profile.SLW.at.ski, profile.MIXPHA.at.ski)
wst.df                    <- data.frame(profile.FZDZ.at.wst, profile.SLW.at.wst, profile.MIXPHA.at.wst)
colnames(sta.df)          <- c("FZDZ1", "FZDZ2", "FZDZ3", "FZDZ4", "FZDZ5", "FZDZ6", "FZDZ7", "FZDZ8", "FZDZ9", "SLW1", "SLW2", "SLW3", "SLW4", "SLW5", "SLW6", "SLW7", "SLW8", "SLW9", "MPHA1", "MPHA2", "MPHA3", "MPHA4", "MPHA5", "MPHA6", "MPHA7", "MPHA8", "MPHA9")
colnames(ski.df)          <- c("FZDZ1", "FZDZ2", "FZDZ3", "FZDZ4", "FZDZ5", "FZDZ6", "FZDZ7", "FZDZ8", "FZDZ9", "SLW1", "SLW2", "SLW3", "SLW4", "SLW5", "SLW6", "SLW7", "SLW8", "SLW9", "MPHA1", "MPHA2", "MPHA3", "MPHA4", "MPHA5", "MPHA6", "MPHA7", "MPHA8", "MPHA9")
colnames(wst.df)          <- c("FZDZ1", "FZDZ2", "FZDZ3", "FZDZ4", "FZDZ5", "FZDZ6", "FZDZ7", "FZDZ8", "FZDZ9", "SLW1", "SLW2", "SLW3", "SLW4", "SLW5", "SLW6", "SLW7", "SLW8", "SLW9", "MPHA1", "MPHA2", "MPHA3", "MPHA4", "MPHA5", "MPHA6", "MPHA7", "MPHA8", "MPHA9")

# TO EXECUTE FOR 3 RadIA ALGOS AND 3 WSMR LOCATIONS:
# 1. SEARCH {ALGO-N}/REPLACE {ALGO-N+1}
# 2. EXECUTE LOOP
# 3. REPEAT FOR N ALGOS
# 4. SWITCH FOR COMMAND BETWEEN 1:16 and 17:31
# 5. REPEAT FOR N ALGOS

# 1 row with 16 columns of plots, 1 for each radar volume
par(mfrow=c(1, 16))
# define margin spacing, starting from x-axis side and moving clockwise
par(mar=c(3, 1, 3, 1))

# user comments/uncomments two lines below to plot first 16 or last 16 times
#for (m in 1:16)  {
for (m in 17:31) {
  
  if (m == 9 | m == 24) {
    plot(profile.MIXPHA.at.wst[m,], NISTkmTOmile(prof.beam.height.wst)*mileTOkft, pch=21, type="b", color="grey", xlim=c(0,1), ylim=c(0.75, 25), xaxt="n", main="MIXPHA")
  } else {
    plot(profile.MIXPHA.at.wst[m,], NISTkmTOmile(prof.beam.height.wst)*mileTOkft, pch=21, type="b", color="grey", xlim=c(0,1), ylim=c(0.75, 25), xaxt="n", yaxt="n", ylab="")
  }  # end of if (m == 9 | m == 24)
  title(paste(substring(radia.SRPC.nc.FRZDZ.filelist.str[m], 1, 2), ":", substring(radia.SRPC.nc.FRZDZ.filelist.str[m], 3, 4), sep=""), line=-67)
  grid()
  
  # initialize these arrays
  status.change <- rep(0,       dim(profile.MIXPHA.at.wst)[2])
  cat.value     <- rep("empty", dim(profile.MIXPHA.at.wst)[2])
  color         <- rep(NA,      dim(profile.MIXPHA.at.wst)[2])
  
  for (n in 1:dim(profile.MIXPHA.at.wst)[2]) {
    
    if (is.na(profile.MIXPHA.at.wst[m,n])) {
      cat.value[n]    <- "empty"
      color[n]        <- NA
    } else {
      if (profile.MIXPHA.at.wst[m,n] < 0.30) {
        cat.value[n]     <- "low"
        color[n]         <- "green"
      } else if (profile.MIXPHA.at.wst[m,n] >= 0.30 & profile.MIXPHA.at.wst[m,n] < 0.60) {
        cat.value[n]     <- "med"
        color[n]         <- "orange"
      } else if (profile.MIXPHA.at.wst[m,n] >= 0.60) {
        cat.value[n]     <- "hi"
        color[n]         <- "red"
      }
    }  # end of if (profile.MIXPHA.at.wst[m,n])
    
    if (n == 1 & cat.value[n] != "empty") {
      print("test")
        status.change[n] <- 1
    }
    
    if (n > 1) {
      if (cat.value[n] != cat.value[n-1]) {
        status.change[n] <- 1
      }
    }  # end of if (n == 1)
    
  }    # end of for (n in 1:dim(profile.DZDZ.at.wst[2]))
  
  ind.status.change     <- which(status.change == 1)
  min.ind.status.change <- min(ind.status.change)
  
  # first, fill in black any heights from the station sfc height up to the first height with radar return
  rect(xleft, wst.alt.ft/1000, xright, NISTkmTOmile(prof.betw.beam.height.wst[min.ind.status.change]) * mileTOkft, density=10, angle=45, col="black", border=NULL, lty=1, lwd=1.5)
  
  # second, fill in height with red/orange/green/black based on presence and magnitude of hazard
  if (length(ind.status.change) > 0) {
    for (kk in 1:(sum(status.change)-1)) {
      ybot   <- NISTkmTOmile(prof.betw.beam.height.wst[ind.status.change[kk]])   * mileTOkft
      ytop   <- NISTkmTOmile(prof.betw.beam.height.wst[ind.status.change[kk+1]]) * mileTOkft
      
      rect(xleft, ybot, xright, ytop, density=10, angle=45, col=color[ind.status.change[kk]], border=NULL, lty=1, lwd=1.5)
    }    # end of for (kk in 1:(sum(status.change)-1))
  }
  
  # third, fill in above highest red/orange/green categorical hazard with black (indicating no radar return)
  rect(xleft, ytop, xright, NISTkmTOmile(7)*mileTOkft, density=10, angle=45, col="black", border=NULL, lty=1, lwd=1)
  
  # finally, fill in from station sfc height down to MSL height=0 with a solid brown rectangle (indicating ground)
  rect(xleft, 0, xright, wst.alt.ft/1000, density=-1, col="chocolate4", border=NULL)
  
}      # end of for (m in 1:16)


