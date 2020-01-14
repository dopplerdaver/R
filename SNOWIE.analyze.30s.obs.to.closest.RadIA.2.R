#----------------------------------------------
#
# Name:    SNOWIE.analyze.30s.obs.to.closest.RadIA.2.R
#
# Purpose: 1) load the 30-sec average UWKA data outputted by SNOWIE.process.30s.obs.to.closest.RadIA.1.R
#          2) make plots and stats
#
# Created: 5.15.2019 dserke
#
#----------------------------------------------

# add libraries
library(lubridate)
library(R.matlab)
library(ncdf4)
library(NISTunits)
library(useful)
library(ggplot2)
library(gridExtra)
library(grid)
library(plotly)
library(RColorBrewer)
library(oceanmap)
library(zoo)
library(R.utils)
py <- plot_ly()

require(dplyr)
require(geosphere)

# remove scientific notation in printing
options(scipen=999)

#----------------------------------------------
# prompt user for input
#----------------------------------------------
inp <- readline(prompt="Enter 1 for single date, 2 for campaign summary: ")
print(paste("user input is: ", inp, sep=""))

#----------------------------------------------
# constants and parameters
#----------------------------------------------
# define yyyymmdd character string to work from
if (inp == 1) {
  radia.yyyy                         <- "2017"
  radia.mm                           <- "03"
  radia.dd                           <- "09"
} else {
  radia.yyyy                         <- "2017"
  radia.mm                           <- "all"
  radia.dd                           <- "SNOWIE"
}
radia.date                         <- paste(radia.yyyy, radia.mm, radia.dd, sep = "")
print(radia.date)

range.gate.dist.km                 <- 0.25

#----------------------------------------------
# define data paths
#----------------------------------------------

# define 10/30s obs data file path
matched.30s.obs.path               <- file.path("/d1/serke/projects/case_studies/SNOWIE/data/RadIA_matchto_AC_data/v3.2/")
matched.30s.obs.filename           <- paste(radia.date, "_RadIA_ints_match.wSLW.v3.csv", sep="")
matched.30s.obs.filename

# define NEXRAD site location csv file path
nexrad.site.dataframetxt.dir       <- file.path("/d1/serke/projects/case_studies/SNOWIE/data/RadIA_data/nexrad_site_data/")

# output data
output.file.dir                    <- file.path("/d1/serke/projects/case_studies/SNOWIE/data/RadIA_matchto_AC_data/")

#----------------------------------------------
# load 10 or 30 second RadIA/AC matched second averaged data set into data frame
#----------------------------------------------
# for 30s obs data, version 2 (w/ koro values)
rawData2                           <- read.csv(paste(matched.30s.obs.path, matched.30s.obs.filename, sep=""), header = TRUE, sep = ",", dec = ".")

# indexing
# for general
#ind.lwc.gt00                           <- which(pmax(rawData2$LWC.nev, rawData2$LWC.ros, na.rm=TRUE) > 0)
ind.d99.gt00                           <- which(rawData2$D99.um > 0)
ind.dmax.gt00                          <- which(rawData2$Dmax.um > 0)
ind.dmax.gt00.cdp.gt01                 <- which(rawData2$Dmax.um > 0 & rawData2$LWC.cdp > 0.01)
ind.LWC.gt01                           <- which(rawData2$LWC.cdp > 0.01)
ind.azim.lt20                          <- which(rawData2$round.azim.radar.to.ac.deg..digits...2. < 20 | rawData2$round.azim.radar.to.ac.deg..digits...2. > 270)
ind.d99.gt00.azim.lt20                 <- intersect(ind.d99.gt00, ind.azim.lt20)
ind.dmax.gt00.azim.lt20                <- intersect(ind.dmax.gt00, ind.azim.lt20)

ind.hi.lwcdiff                         <- which(rawData2$LWC.ros - rawData2$LWC.cdp >= 0.03)   
median.hi.lwcdiff                      <- median(rawData2$FZDZ.int.match.matrix.max[ind.hi.lwcdiff], na.rm=TRUE)
ind.lo.lwc.diff                        <- which(rawData2$LWC.ros - rawData2$LWC.cdp < 0.04)

#---------------------
# all stuff to next dashed line was coded after talking w sarah on 6/24/2019
#ind.cdp.gt01                            <- which(rawData2$LWC.cdp > 0.01)

thresh.matrix                           <- seq(0, 1.00, 0.02)

#   for FZDZ
ind.dmax.gt100.cdp.gt01.mlf.gt20        <- which(rawData2$Dmax.um >= 100 & rawData2$LWC.cdp > 0.01 & rawData2$mass.liq.frac > 0.20)
# choose the hard criteria for FZDZ....
ind.ac.appO                             <- intersect(ind.hi.lwcdiff, ind.dmax.gt100.cdp.gt01.mlf.gt20)
# or choose the soft criteria for FZDZ ...
#ind.ac.appO                             <- ind.dmax.gt100.cdp.gt01.mlf.gt20
ind.dmax.lt100.cdp.gt01.mlf.gt20        <- which(rawData2$Dmax.um < 100 & rawData2$LWC.cdp > 0.01 & rawData2$mass.liq.frac > 0.20)
# choose the hard criteria for small-drop ....
ind.ac.appC                             <- intersect(ind.lo.lwc.diff, ind.dmax.lt100.cdp.gt01.mlf.gt20)
# or choose the soft criteria for FZDZ ...
#ind.ac.appC                             <- ind.dmax.lt100.cdp.gt01.mlf.gt20
PODY.FZDZ                               <- array(0, length(thresh.matrix))
PODN.FZDZ                               <- array(0, length(thresh.matrix))
FAR.FZDZ                                <- array(0, length(thresh.matrix))
MISS.FZDZ                               <- array(0, length(thresh.matrix))
SPEC.FZDZ                               <- array(0, length(thresh.matrix))
for (n in 1:length(thresh.matrix)) {
# for (n in 60:60) {
  print("-------")
  print(thresh.matrix[n])
  #ind.dmax.gt100.cdp.gt01.mlf.gt80        <- which(rawData2$Dmax.um >= 100 & rawData2$LWC.cdp > 0.01 & rawData2$mass.liq.frac > 0.80)
  #ind.dmax.gt100.cdp.gt01.mlf.gt80.FZ.NA  <- which(rawData2$Dmax.um >= 100 & rawData2$LWC.cdp > 0.01 & rawData2$mass.liq.frac > 0.80 & is.na(rawData2$FZDZ.int.match.matrix.max))
  num.dmax.gt100.cdp.gt01.mlf.gt20        <- length(ind.dmax.gt100.cdp.gt01.mlf.gt20)
  num.ac.appO                             <- length(ind.ac.appO)                                                                 # denom of FAR for SLW, denom of PODN for SLW
  num.ac.appO.notNA                       <- length(which(!is.na(rawData2$FZDZ.int.match.matrix.max[ind.ac.appO])))              # denom of PODY for FZDZ, denom for MISSR for FZDZ
  num.ac.appO.FZDZgtval                   <- length(which(rawData2$FZDZ.int.match.matrix.max[ind.ac.appO] >= thresh.matrix[n]))  # numer of PODY for FZDZ
  num.ac.appC.FZDZgtval                   <- length(which(rawData2$FZDZ.int.match.matrix.max[ind.ac.appC] >= thresh.matrix[n]))  # numer of FAR for FZDZ
  num.ac.appO.FZDZltval                   <- length(which(rawData2$FZDZ.int.match.matrix.max[ind.ac.appC] <  thresh.matrix[n]))  # numer of MISS for FZDZ
  num.ac.appC.FZDZltval                   <- length(which(rawData2$FZDZ.int.match.matrix.max[ind.ac.appO] <  thresh.matrix[n]))  # numer of PODN for FZDZ
  PODY.FZDZ[n]                            <- num.ac.appO.FZDZgtval / num.ac.appO.notNA
  PODN.FZDZ[n]                            <- num.ac.appC.FZDZltval / num.ac.appO
  FAR.FZDZ[n]                             <- num.ac.appC.FZDZgtval / num.ac.appO
  MISS.FZDZ[n]                            <- num.ac.appO.FZDZltval / num.ac.appO.notNA
  SPEC.FZDZ[n]                            <- 1 - FAR.FZDZ[n]
  #num.dmax.gt100.cdp.gt01.mlf.gt20
  print(num.ac.appO)
  print(num.ac.appO.notNA)
  print(num.ac.appO.FZDZgtval)
  #print(num.ac.appC.FZDZgtval)
  #print(num.ac.appO.FZDZltval)
  #print(num.ac.appC.FZDZltval)
  print(PODY.FZDZ[n])
  print(PODN.FZDZ[n])
  print(FAR.FZDZ[n])
  print(MISS.FZDZ[n])
  print(SPEC.FZDZ[n])
}
PODY.FZDZ.extend <- c(1, PODY.FZDZ, 0)
FAR.FZDZ.extend  <- c(1, FAR.FZDZ, 0)

#   for SLW
ind.dmax.lt100.cdp.gt01.mlf.gt20        <- which(rawData2$Dmax.um < 100 & rawData2$LWC.cdp > 0.01 & rawData2$mass.liq.frac > 0.20)
ind.ac.appC                             <- intersect(ind.lo.lwc.diff, ind.dmax.lt100.cdp.gt01.mlf.gt20) 
PODY.SLW                                <- array(0, length(thresh.matrix))
PODN.SLW                                <- array(0, length(thresh.matrix))
FAR.SLW                                 <- array(0, length(thresh.matrix))
MISS.SLW                                <- array(0, length(thresh.matrix))
SPEC.SLW                                <- array(0, length(thresh.matrix))
for (m in 1:length(thresh.matrix)) {
# for (m in 60:60) {
  num.ac.appC                             <- length(ind.ac.appC)                                                                # denom of FAR for FZDZ, denom of PODN for FZDZ
  #num.ac.appC
  num.ac.appC.notNA                       <- length(which(!is.na(rawData2$SLW.int.match.column.max[ind.ac.appC])))              # denom of PODY for SLW, denom for MISSR for SLW, denom of PODN for FZDZ
  #num.ac.appC.notNA
  num.ac.appC.SLWgtval                    <- length(which(rawData2$SLW.int.match.column.max[ind.ac.appC] >= thresh.matrix[m]))  # numer of PODY for SLW
  #num.ac.appC.SLWgtval
  num.ac.appO.SLWgtval                    <- length(which(rawData2$SLW.int.match.column.max[ind.ac.appO] >= thresh.matrix[m]))  # numer of FAR for SLW
  #num.ac.appO.SLWgtval
  num.ac.appC.SLWltval                    <- length(which(rawData2$SLW.int.match.column.max[ind.ac.appO] <  thresh.matrix[m]))  # numer of MISS for SLW
  #num.ac.appC.SLWltval
  num.ac.appO.SLWltval                    <- length(which(rawData2$SLW.int.match.column.max[ind.ac.appC] <  thresh.matrix[m]))  # numer of PODN for SLW
  #num.ac.appO.SLWltval

  num.dmax.lt100.cdp.gt01.mlf.gt20        <- length(ind.dmax.lt100.cdp.gt01.mlf.gt20)
  #num.dmax.lt100.cdp.gt01.mlf.gt20
  
  PODY.SLW[m]                            <- num.ac.appC.SLWgtval / num.ac.appC.notNA
  print(PODY.SLW[m])
  PODN.SLW[m]                            <- num.ac.appO.SLWltval / num.ac.appC
  print(PODN.SLW[m])
  FAR.SLW[m]                             <- num.ac.appO.SLWgtval / num.ac.appC
  print(FAR.SLW[m])
  MISS.SLW[m]                            <- num.ac.appC.SLWltval / num.ac.appC.notNA
  print(MISS.SLW[m])
  SPEC.SLW[m]                            <- 1 - FAR.SLW[m]
  print(SPEC.SLW[m])
}
PODY.SLW.extend <- c(1, PODY.SLW, 0)
FAR.SLW.extend  <- c(1, FAR.SLW, 0)

#   for MIXPHA
# choose one of the following two lines to uncomment
# for AMS RADAR poster, ran next line for LWC-CDP 0., 0.15, and 0.25
ind.ac.mixpha                           <- which(rawData2$LWC.cdp >= 0.05 & rawData2$mass.liq.frac > 0.20 & rawData2$mass.liq.frac < 0.80 & rawData2$MIXPHA.int.match.matrix.max > 0.33)
#ind.ac.mixpha                           <- which(rawData2$LWC.cdp > 0.01 & rawData2$mass.liq.frac > 0.20 & rawData2$mass.liq.frac < 0.80)
#ind.ac.mixpha                           <- which(rawData2$LWC.cdp > 0.01 & rawData2$mass.liq.frac > 0.20 & rawData2$mass.liq.frac < 0.80 | rawData2$mass.liq.frac == 0 | rawData2$mass.liq.frac == 1.00)
print(length(ind.ac.mixpha))
ind.ac.cdp.gt01.mlf.gt95                <- which(rawData2$LWC.cdp >= 0.01 & rawData2$mass.liq.frac > 0.95 & rawData2$Dmax.um > 0)
ind.ac.notmixpha                        <- union(ind.ac.iceonly3, ind.ac.cdp.gt01.mlf.gt95)
PODY.MPHA                                <- array(0, length(thresh.matrix))
PODN.MPHA                                <- array(0, length(thresh.matrix))
FAR.MPHA                                 <- array(0, length(thresh.matrix))
MISS.MPHA                                <- array(0, length(thresh.matrix))
SPEC.MPHA                                <- array(0, length(thresh.matrix))
for (p in 1:length(thresh.matrix)) {
# for (p in 60:60) {
  num.ac.MPHA                             <- length(ind.ac.mixpha)                                                      
  print(num.ac.MPHA)
  num.ac.notMPHA                          <- length(ind.ac.notmixpha)                                                                   # denom of FAR for MIXPHA, denom for PODN for MIXPHA
  print(num.ac.notMPHA)
  num.ac.MPHA.MPHA.gtval                  <- length(which(rawData2$MIXPHA.int.match.column.max[ind.ac.mixpha]    >  thresh.matrix[p]))  # denom of PODY for MIXPHA, denom of MISSR for MIXPHA
  print(num.ac.MPHA.MPHA.gtval)
  num.ac.MPHA.MPHA.gtval                  <- length(which(rawData2$MIXPHA.int.match.column.max[ind.ac.mixpha]    >= thresh.matrix[p]))  # numer of PODY for MIXPHA
  #num.ac.MPHA.MPHA.gtval
  num.ac.MPHA.MPHA.ltval                  <- length(which(rawData2$MIXPHA.int.match.column.max[ind.ac.mixpha]     < thresh.matrix[p] & rawData2$MIXPHA.int.match.column.max[ind.ac.mixpha] > 0.33))  # numer of MISSR for MIXPHA
  #num.ac.MPHA.MPHA.ltval
  num.ac.notMPHA.MPHA.gtval               <- length(which(rawData2$MIXPHA.int.match.column.max[ind.ac.notmixpha] >= thresh.matrix[p]))  # numer of FAR for MIXPHA
  #num.ac.notMPHA.MPHA.gtval
  num.ac.notMPHA.MPHA.ltval               <- length(which(rawData2$MIXPHA.int.match.column.max[ind.ac.notmixpha] <  thresh.matrix[p]))  # numer of PODN for MIXPHA
  #num.ac.notMPHA.MPHA.ltval
  
  PODY.MPHA[p]                            <- num.ac.MPHA.MPHA.gtval    / num.ac.MPHA
  #print(PODY.MPHA[p])
  PODN.MPHA[p]                            <- num.ac.notMPHA.MPHA.ltval / num.ac.notMPHA
  #print(PODN.MPHA[p])
  FAR.MPHA[p]                             <- num.ac.notMPHA.MPHA.gtval / num.ac.notMPHA
  #print(FAR.MPHA[p])
  MISS.MPHA[p]                            <- num.ac.MPHA.MPHA.ltval    / num.ac.MPHA.MPHA.gtval
  #print(MISS.MPHA[p])
  SPEC.MPHA[p]                            <- 1 - FAR.MPHA[n]
  #print(SPEC.MPHA[p])
} 
     
#   for ICEONLY
thresh.matrix.iceonly                    <- seq(0.30, 1.00, 0.02)
#ind.ac.iceonly                          <- which(rawData2$LWC.cdp <= 0 & rawData2$LWC.nev <= 0 & rawData2$LWC.ros <= 0 & rawData2$mass.liq.frac <= 0.05)
#ind.ac.iceonly.keep                     <- c(4, 5, 6, 7, 20, 21, 22, 23, 24, 25, 26, 76, 77, 78, 79, 80, 81, 82, 89, 90, 91, 92, 93, 94, 95, 96, 97, 98, 112, 113, 114, 115, 116, 117, 118, 119, 120, 121, 122, 123, 124, 125, 126, 134, 135, 136, 150, 151, 152, 153, 154, 155, 156, 157, 158, 159, 190, 191, 192, 193, 194, 195, 196, 197, 198, 199, 200, 201, 202, 203, 204, 205, 206, 207, 208, 209, 210, 211, 212, 213, 227, 228, 229, 230, 231, 232, 273, 274, 275, 276, 277, 278)  # indices with at least a minute of ice-only data on either side of the given index
#length(ind.ac.iceonly[ind.ac.iceonly.keep])
#ind.ac.iceonly2                         <- which(rawData2$LWC.cdp <= 0 & rawData2$LWC.nev <= 0 & rawData2$LWC.ros <= 0 & rawData2$mass.liq.frac <= 0.00 & rawData2$SLW.int.match.column.max <= 0.75)
#length(ind.ac.iceonly)
ind.ac.iceonly3                         <- which(rawData2$LWC.cdp == 0 & rawData2$Dmax.um == 0 & !is.na(rawData2$FZDZ.int.match.column.max) & !is.na(rawData2$SLW.int.match.column.max) & !is.na(rawData2$MIXPHA.int.match.column.max))
print(length(ind.ac.iceonly3))
ind.ac.iceonlydetect                    <- which(rawData2$LWC.cdp == 0 & rawData2$Dmax.um == 0 & rawData2$FZDZ.int.match.column.max <= 0.60 & rawData2$SLW.int.match.column.max < 0.75 & rawData2$MIXPHA.int.match.column.max <= 0.60)
print(length(ind.ac.iceonlydetect))
ind.ac.noticeonlydetect                 <- which(rawData2$LWC.cdp > 0 & rawData2$Dmax.um > 0 & rawData2$FZDZ.int.match.column.max <= 0.60 & rawData2$SLW.int.match.column.max < 0.75 & rawData2$MIXPHA.int.match.column.max <= 0.60)
print(length(ind.ac.noticeonlydetect))
ind.ac.noticeonlydetect.notNA           <- which(rawData2$LWC.cdp > 0 & rawData2$Dmax.um > 0 & !is.na(rawData2$FZDZ.int.match.column.max) & !is.na(rawData2$SLW.int.match.column.max) & !is.na(rawData2$MIXPHA.int.match.column.max))
print(length(ind.ac.noticeonlydetect.notNA))
ind.ac.noticeonlydetect.noticeonly      <- which(rawData2$LWC.cdp > 0 & rawData2$Dmax.um > 0 & rawData2$FZDZ.int.match.column.max > 0.60 & rawData2$SLW.int.match.column.max >= 0.75 & rawData2$MIXPHA.int.match.column.max > 0.60)
print(length(ind.ac.noticeonlydetect.noticeonly))
ind.all                                 <- which(rawData2$LWC.cdp == rawData2$LWC.cdp)
ind.ac.noticeonly3                      <- setdiff(ind.all, ind.ac.iceonly3)                      # denom of FAR for ALLICE, denom of PODN for ALLICE
length(ind.ac.noticeonly3)
#num.ac.noticeonly.ALLICE                <-
test1000                                <- which(rawData2$LWC.cdp > 0 & rawData2$Dmax > 0)
length(test1000)
ind.ac.noticeonly.notALLICE             <- which(rawData2$LWC.cdp > 0 & rawData2$Dmax > 0 & (rawData2$FZDZ.int.match.column.max > 0.60 | rawData2$SLW.int.match.column.max > 0.75 | rawData2$MIXPHA.int.match.column.max > 0.62))
length(ind.ac.noticeonly.notALLICE)
PODY.ICEONLY                            <- array(0, length(thresh.matrix))
PODN.ICEONLY                            <- array(0, length(thresh.matrix))
FAR.ICEONLY                             <- array(0, length(thresh.matrix))
MISS.ICEONLY                            <- array(0, length(thresh.matrix))
SPEC.ICEONLY                            <- array(0, length(thresh.matrix))
for (q in 1:length(thresh.matrix)) {
  # for (p in 60:60) {
  num.ac.ICEONLY                          <- length(ind.ac.iceonlydetect)                                                      
  print(num.ac.ICEONLY)
  num.ac.notICEONLY                       <- length(ind.ac.noticeonlydetect)                                                                   # denom of FAR for MIXPHA, denom for PODN for MIXPHA
  print(num.ac.notICEONLY)
  num.ac.ICEONLY.ICEONLY.gtval            <- length(which(rawData2$FZDZ.int.match.column.max[ind.ac.iceonlydetect]    >  thresh.matrix[q]-0.46 & rawData2$SLW.int.match.column.max[ind.ac.iceonlydetect]    >  thresh.matrix[q] & rawData2$MIXPHA.int.match.column.max[ind.ac.iceonlydetect] >  thresh.matrix[q]-0.34))  # denom of PODY for MIXPHA, denom of MISSR for MIXPHA
  print(num.ac.ICEONLY.ICEONLY.gtval)
  num.ac.ICEONLY.ICEONLY.gtval            <- length(which(rawData2$FZDZ.int.match.column.max[ind.ac.iceonlydetect]    >= thresh.matrix[q]-0.46 & rawData2$SLW.int.match.column.max[ind.ac.iceonlydetect]    >  thresh.matrix[q] & rawData2$MIXPHA.int.match.column.max[ind.ac.iceonlydetect] >  thresh.matrix[q]-0.34))  # numer of PODY for MIXPHA
  #num.ac.ICEONLY.ICEONLY.gtval
  num.ac.ICEONLY.ICEONLY.ltval            <- length(which(rawData2$FZDZ.int.match.column.max[ind.ac.iceonlydetect]     < thresh.matrix[q]-0.46 & rawData2$SLW.int.match.column.max[ind.ac.iceonlydetect]    >  thresh.matrix[q] & rawData2$MIXPHA.int.match.column.max[ind.ac.iceonlydetect] >  thresh.matrix[q]-0.34 & rawData2$MIXPHA.int.match.column.max[ind.ac.iceonlydetect] > 0.33))  # numer of MISSR for MIXPHA
  #num.ac.ICEONLY.ICEONLY.ltval
  num.ac.notICEONLY.ICEONLY.gtval         <- length(which(rawData2$FZDZ.int.match.column.max[ind.ac.noticeonlydetect] >= thresh.matrix[q]-0.46 & rawData2$SLW.int.match.column.max[ind.ac.iceonlydetect]    >  thresh.matrix[q] & rawData2$MIXPHA.int.match.column.max[ind.ac.iceonlydetect] >  thresh.matrix[q]-0.34))  # numer of FAR for MIXPHA
  #num.ac.notICEONLY.ICEONLY.gtval
  num.ac.notICEONLY.ICEONLY.ltval         <- length(which(rawData2$FZDZ.int.match.column.max[ind.ac.noticeonlydetect] <  thresh.matrix[q]-0.46 & rawData2$SLW.int.match.column.max[ind.ac.iceonlydetect]    >  thresh.matrix[q] & rawData2$MIXPHA.int.match.column.max[ind.ac.iceonlydetect] >  thresh.matrix[q]-0.34))  # numer of PODN for MIXPHA
  #num.ac.notICEONLY.ICEONLY.ltval
  
  PODY.ICEONLY[q]                         <- num.ac.ICEONLY.ICEONLY.gtval / num.ac.ICEONLY
  #print(PODY.ICEONLY[q])
  PODN.ICEONLY[q]                         <- num.ac.notMPHA.ICEONLY.ltval / num.ac.notICEONLY
  #print(PODN.ICEONLY[q])
  FAR.ICEONLY[q]                          <- num.ac.notMPHA.ICEONLY.gtval / num.ac.notICEONLY
  #print(FAR.ICEONLY[q])
  MISS.ICEONLY[q]                         <- num.ac.ICEONLY.ICEONLY.ltval / num.ac.ICEONLY.ICEONLY.gtval
  #print(MISS.ICEONLY[q])
  SPEC.ICEONLY[q]                         <- 1 - FAR.ICEONLY[n]
  #print(SPEC.ICEONLY[q])
} 
                                           
#--------------------------

FZDZ.mean.ind.ac.appO   <- mean(rawData2$FZDZ.int.match.matrix.max[ind.ac.appO], na.rm=TRUE)
FZDZ.median.ind.ac.appO <- median(rawData2$FZDZ.int.match.matrix.max[ind.ac.appO], na.rm=TRUE)

ind.dmax.gt00.FZDZ.gt.60                <- which(rawData2$Dmax.um >  0 & rawData2$Dmax.um <= 500 & rawData2$FZDZ.int.match.matrix.max >= 0.60)
ind.dmax.gt00.FZDZ.gt.00                <- which(rawData2$Dmax.um >  0 & rawData2$Dmax.um <= 500 & rawData2$FZDZ.int.match.matrix.max >= 0.00)
#ind.dmax.gt00.FZDZ.gt.60.lwcdiff.gt.03  <- which(rawData2$dmax.um >  0 & rawData2$dmax.um <= 500 & rawData2$FZDZ.int.match.matrix.max >= 0.60 & (rawData2$LWC.ros - rawData2$LWC.nev) > 0.05)
#ind.dmax.gt00.FZDZ.gt.00.lwcdiff.gt.03  <- which(rawData2$dmax.um >  0 & rawData2$dmax.um <= 500 & rawData2$FZDZ.int.match.matrix.max >= 0.00 & (rawData2$LWC.ros - rawData2$LWC.nev) > 0.05)
#ind.dmax.gt00.FZDZ.gt.60.lwcdiff.lt.03  <- which(rawData2$dmax.um >  0 & rawData2$dmax.um <= 500 & rawData2$FZDZ.int.match.matrix.max >= 0.60 & (rawData2$LWC.ros - rawData2$LWC.nev) <= 0.05)
#ind.dmax.gt00.FZDZ.gt.00.lwcdiff.lt.03  <- which(rawData2$dmax.um >  0 & rawData2$dmax.um <= 500 & rawData2$FZDZ.int.match.matrix.max >= 0.00 & (rawData2$LWC.ros - rawData2$LWC.nev) <= 0.05)
ind.dmax.gt00.FZDZ.eq.NA                <- which(rawData2$dmax.um >  0 & rawData2$dmax.um <= 500 & is.na(rawData2$FZDZ.int.match.matrix.max))
#ind.dmax.gt00.FZDZ.eq.NA.lwcdiff.gt.03  <- which(rawData2$dmax.um >  0 & rawData2$dmax.um <= 500 & is.na(rawData2$FZDZ.int.match.matrix.max)  & (rawData2$LWC.ros - rawData2$LWC.nev) > 0.05)
#ind.dmax.gt00.FZDZ.eq.NA.lwcdiff.lt.03  <- which(rawData2$dmax.um >  0 & rawData2$dmax.um <= 500 & is.na(rawData2$FZDZ.int.match.matrix.max)  & (rawData2$LWC.ros - rawData2$LWC.nev) <= 0.05)
ind.dmax.eq00.FZDZ.eq.NA                <- which(rawData2$dmax.um == 0 & rawData2$dmax.um <= 500 & is.na(rawData2$FZDZ.int.match.matrix.max))
ind.dmax.eq00.FZDZ.gt.60                <- which(rawData2$dmax.um == 0 & rawData2$dmax.um <= 500 & rawData2$FZDZ.int.match.matrix.max >= 0.60)
ind.dmax.eq00.FZDZ.gt.00                <- which(rawData2$dmax.um == 0 & rawData2$dmax.um <= 500 & rawData2$FZDZ.int.match.matrix.max > 0)
ind.FZDZ.gt.00                         <- intersect(ind.azim.lt20, which(rawData2$FZDZ.int.match.matrix.max >  0.00))
ind.FZDZ.gt.60                         <- intersect(ind.azim.lt20, which(rawData2$FZDZ.int.match.matrix.max >= 0.60))
ind.dmax.eq00                           <- intersect(ind.azim.lt20, which(rawData2$Dmax.um == 0))
ind.dmax.eq00.FZDZ.eqNA                 <- intersect(ind.azim.lt20, which(rawData2$Dmax.um == 0 & is.na(rawData2$FZDZ.int.match.matrix.max)))
ind.dmax.gt00                           <- intersect(ind.azim.lt20, which(rawData2$Dmax.um > 0))
ind.dmax.gt50                           <- intersect(ind.azim.lt20, which(rawData2$Dmax.um > 50 & rawData2$dmax.um <= 500))
ind.dmax.gt50.FZDZ.eq.NA                <- intersect(ind.azim.lt20, which(rawData2$Dmax.um > 50 & rawData2$dmax.um <= 500 & is.na(rawData2$FZDZ.int.match.matrix.max)))
ind.FZDZ.eq.NA                         <- intersect(ind.azim.lt20, which(is.na(rawData2$FZDZ.int.match.matrix.max)))

#   for MPHA
ind.dmax.gt00.MPHA.gt.60                <- which(rawData2$Dmax.um >  0 & rawData2$Dmax.um <= 500 & rawData2$MIXPHA.int.match.matrix.max >= 0.60)
ind.dmax.gt00.MPHA.gt.00                <- which(rawData2$Dmax.um >  0 & rawData2$Dmax.um <= 500 & rawData2$MIXPHA.int.match.matrix.max > 0.00)
ind.dmax.gt00.MPHA.gt.60.lwcdiff.gt.03  <- which(rawData2$Dmax.um >  0 & rawData2$Dmax.um <= 500 & rawData2$MIXPHA.int.match.matrix.max >= 0.60 & (rawData2$LWC.ros - rawData2$LWC.nev) > 0.05)
ind.dmax.gt00.MPHA.gt.00.lwcdiff.gt.03  <- which(rawData2$Dmax.um >  0 & rawData2$Dmax.um <= 500 & rawData2$MIXPHA.int.match.matrix.max > 0.00  & (rawData2$LWC.ros - rawData2$LWC.nev) > 0.05)
ind.dmax.gt00.MPHA.gt.60.lwcdiff.lt.03  <- which(rawData2$Dmax.um >  0 & rawData2$Dmax.um <= 500 & rawData2$MIXPHA.int.match.matrix.max >= 0.60 & (rawData2$LWC.ros - rawData2$LWC.nev) <= 0.05)
ind.dmax.gt00.MPHA.gt.00.lwcdiff.lt.03  <- which(rawData2$Dmax.um >  0 & rawData2$Dmax.um <= 500 & rawData2$MIXPHA.int.match.matrix.max > 0.00  & (rawData2$LWC.ros - rawData2$LWC.nev) <= 0.05)
ind.dmax.gt00.MPHA.eq.NA                <- which(rawData2$Dmax.um >  0 & rawData2$Dmax.um <= 500 & is.na(rawData2$MIXPHA.int.match.matrix.max))
ind.dmax.eq00.MPHA.eq.NA                <- which(rawData2$Dmax.um == 0 & rawData2$Dmax.um <= 500 & is.na(rawData2$MIXPHA.int.match.matrix.max))
ind.dmax.eq00.MPHA.gt.60                <- which(rawData2$Dmax.um == 0 & rawData2$Dmax.um <= 500 & rawData2$MIXPHA.int.match.matrix.max >= 0.60)
ind.dmax.eq00.MPHA.gt.00                <- which(rawData2$Dmax.um == 0 & rawData2$Dmax.um <= 500 & rawData2$MIXPHA.int.match.matrix.max > 0)

#   for SLW
ind.dmax.gt00.SLW.gt.60                <- which(rawData2$Dmax.um >  0 & rawData2$Dmax.um <= 500 & rawData2$SLW.int.match.matrix.max >= 0.60)
ind.dmax.gt00.SLW.gt.00                <- which(rawData2$Dmax.um >  0 & rawData2$Dmax.um <= 500 & rawData2$SLW.int.match.matrix.max > 0.00)
ind.dmax.gt00.SLW.gt.60.lwcdiff.gt.03  <- which(rawData2$Dmax.um >  0 & rawData2$Dmax.um <= 500 & rawData2$SLW.int.match.matrix.max >= 0.60 & (rawData2$LWC.ros - rawData2$LWC.nev) > 0.05)
ind.dmax.gt00.SLW.gt.00.lwcdiff.gt.03  <- which(rawData2$Dmax.um >  0 & rawData2$Dmax.um <= 500 & rawData2$SLW.int.match.matrix.max > 0.00  & (rawData2$LWC.ros - rawData2$LWC.nev) > 0.05)
ind.dmax.gt00.SLW.gt.60.lwcdiff.lt.03  <- which(rawData2$Dmax.um >  0 & rawData2$Dmax.um <= 500 & rawData2$SLW.int.match.matrix.max >= 0.60 & (rawData2$LWC.ros - rawData2$LWC.nev) <= 0.05)
ind.dmax.gt00.SLW.gt.00.lwcdiff.lt.03  <- which(rawData2$Dmax.um >  0 & rawData2$Dmax.um <= 500 & rawData2$SLW.int.match.matrix.max > 0.00  & (rawData2$LWC.ros - rawData2$LWC.nev) <= 0.05)
ind.dmax.gt00.SLW.eq.NA                <- which(rawData2$Dmax.um >  0 & rawData2$Dmax.um <= 500 & is.na(rawData2$SLW.int.match.matrix.max))
ind.dmax.eq00.SLW.eq.NA                <- which(rawData2$Dmax.um == 0 & rawData2$Dmax.um <= 500 & is.na(rawData2$SLW.int.match.matrix.max))
ind.dmax.eq00.SLW.gt.60                <- which(rawData2$Dmax.um == 0 & rawData2$Dmax.um <= 500 & rawData2$SLW.int.match.matrix.max >= 0.60)
ind.dmax.eq00.SLW.gt.00                <- which(rawData2$Dmax.um == 0 & rawData2$Dmax.um <= 500 & rawData2$SLW.int.match.matrix.max > 0)

#--------------------
# length of indices

# for general
num.total.pts                         <- dim(rawData2)[1]
#num.lwc.gt00                          <- length(ind.lwc.gt00)
num.LWC.gt01                          <- length(ind.LWC.gt01)
num.dmax.gt00                          <- length(ind.dmax.gt00)
num.dmax.gt00.cdp.gt01                <- length(ind.dmax.gt00.cdp.gt01)

#   for FZDZ
#num.FZDZ.gt00.dmax.gt00.cdp.gt01       <- length(ind.FZDZ.gt00.dmax.gt00.cdp.gt01)
#num.FZDZ.gt60.dmax.gt00.cdp.gt01       <- length(ind.FZDZ.gt60.dmax.gt00.cdp.gt01)
num.dmax.gt00.FZDZ.gt.60               <- length(ind.dmax.gt00.FZDZ.gt.60)
num.dmax.gt00.FZDZ.gt.00               <- length(ind.dmax.gt00.FZDZ.gt.00)
#num.dmax.gt00.FZDZ.gt.60.lwcdiff.gt.03 <- length(ind.dmax.gt00.FZDZ.gt.60.lwcdiff.gt.03)
#num.dmax.gt00.FZDZ.gt.00.lwcdiff.gt.03 <- length(ind.dmax.gt00.FZDZ.gt.00.lwcdiff.gt.03)
#num.dmax.gt00.FZDZ.gt.60.lwcdiff.lt.03 <- length(ind.dmax.gt00.FZDZ.gt.60.lwcdiff.lt.03)
#num.dmax.gt00.FZDZ.gt.00.lwcdiff.lt.03 <- length(ind.dmax.gt00.FZDZ.gt.00.lwcdiff.lt.03)
num.dmax.gt00.FZDZ.eq.NA               <- length(ind.dmax.gt00.FZDZ.eq.NA)
#num.dmax.gt00.FZDZ.eq.NA.lwcdiff.gt.03 <- length(ind.dmax.gt00.FZDZ.eq.NA.lwcdiff.gt.03)
#num.dmax.gt00.FZDZ.eq.NA.lwcdiff.lt.03 <- length(ind.dmax.gt00.FZDZ.eq.NA.lwcdiff.lt.03)
#num.dmax.eq00.FZDZ.gt.60               <- length(ind.dmax.eq00.FZDZ.gt.60)
num.dmax.eq00.FZDZ.gt.00               <- length(ind.dmax.eq00.FZDZ.gt.00)
num.dmax.eq00.FZDZ.eq.NA               <- length(ind.dmax.eq00.FZDZ.eq.NA)

num.dmax.gt00.FZDZ.gt.60
num.dmax.gt00.FZDZ.gt.00
#num.dmax.gt00.FZDZ.gt.60.lwcdiff.gt.03
#num.dmax.gt00.FZDZ.gt.00.lwcdiff.gt.03
#num.dmax.gt00.FZDZ.gt.60.lwcdiff.lt.03
#num.dmax.gt00.FZDZ.gt.00.lwcdiff.lt.03
num.dmax.gt00.FZDZ.eq.NA 
#num.dmax.gt00.FZDZ.eq.NA.lwcdiff.gt.03
#num.dmax.gt00.FZDZ.eq.NA.lwcdiff.lt.03

num.FZDZ.gt.00                        <- length(ind.FZDZ.gt.00)
num.FZDZ.gt.60                        <- length(ind.FZDZ.gt.60)
num.dmax.eq00                          <- length(ind.dmax.eq00)
num.dmax.eq00.FZDZ.eqNA                <- length(ind.dmax.eq00.FZDZ.eqNA)
num.dmax.gt50                          <- length(ind.dmax.gt50)
num.dmax.gt50.FZDZ.eqNA                <- length(ind.dmax.gt50.FZDZ.eq.NA)
num.FZDZ.eqNA                         <- length(ind.FZDZ.eq.NA)

#   for MPHA
num.dmax.gt00.MPHA.gt.60  <- length(ind.dmax.gt00.MPHA.gt.60)
num.dmax.gt00.MPHA.gt.00  <- length(ind.dmax.gt00.MPHA.gt.00)
num.dmax.gt00.MPHA.eq.NA  <- length(ind.dmax.gt00.MPHA.eq.NA)
num.dmax.eq00.MPHA.gt.60  <- length(ind.dmax.eq00.MPHA.gt.60)
num.dmax.eq00.MPHA.gt.00  <- length(ind.dmax.eq00.MPHA.gt.00)
num.dmax.eq00.MPHA.eq.NA  <- length(ind.dmax.eq00.MPHA.eq.NA)

#   for SLW
num.dmax.gt00.SLW.gt.60  <- length(ind.dmax.gt00.SLW.gt.60)
num.dmax.gt00.SLW.gt.00  <- length(ind.dmax.gt00.SLW.gt.00)
num.dmax.gt00.SLW.eq.NA  <- length(ind.dmax.gt00.SLW.eq.NA)
num.dmax.eq00.SLW.gt.60  <- length(ind.dmax.eq00.SLW.gt.60)
num.dmax.eq00.SLW.gt.00  <- length(ind.dmax.eq00.SLW.gt.00)
num.dmax.eq00.SLW.eq.NA  <- length(ind.dmax.eq00.SLW.eq.NA)

# create num cases and num total 30-s points
num.total                <- length(rawData2$Dmax.um)

out                      <- aggregate(data.frame(count=rawData2$radia.date), list(value=rawData2$radia.date), length)
ind.cases                <- which(out$value > 20170000) # value in out df should be greater than this number if the yyyymmdd has been processed
num.cases                <- length(ind.cases)           # num of cases that have been processed

#-------------------------
# create df for dmax vs FZDZ max ints: number of incidences
#-------------------------
#   for FZDZ
#Dmaxgt0.FZDZ              <- c(num.dmax.gt00.FZDZ.gt.60, num.dmax.gt00.FZDZ.gt.00, num.dmax.gt00.FZDZ.eq.NA)
#Dmaxeq0.FZDZ              <- c(num.dmax.eq00.FZDZ.gt.60, num.dmax.eq00.FZDZ.gt.00, num.dmax.eq00.FZDZ.eq.NA)
#Dmax.FZDZ.df              <- data.frame(Dmaxgt0.FZDZ, Dmaxeq0.FZDZ)
#rownames(Dmax.FZDZ.df)    <- c("FZDZgt60", "FZDZgt0", "FZDZeqNA")

#   for MPHA
#Dmaxgt0.MPHA              <- c(num.dmax.gt00.MPHA.gt.60, num.dmax.gt00.MPHA.gt.00, num.dmax.gt00.MPHA.eq.NA)
#Dmaxeq0.MPHA              <- c(num.dmax.eq00.MPHA.gt.60, num.dmax.eq00.MPHA.gt.00, num.dmax.eq00.MPHA.eq.NA)
#Dmax.MPHA.df              <- data.frame(Dmaxgt0.MPHA, Dmaxeq0.MPHA)
#rownames(Dmax.MPHA.df)    <- c("MPHAgt60", "MPHAgt0", "MPHAeqNA")

#   for SLW
#Dmaxgt0.SLW              <- c(num.dmax.gt00.SLW.gt.60, num.dmax.gt00.SLW.gt.00, num.dmax.gt00.SLW.eq.NA)
#Dmaxeq0.SLW              <- c(num.dmax.eq00.SLW.gt.60, num.dmax.eq00.SLW.gt.00, num.dmax.eq00.SLW.eq.NA)
#Dmax.SLW.df              <- data.frame(Dmaxgt0.SLW, Dmaxeq0.SLW)
#rownames(Dmax.SLW.df)    <- c("SLWgt60", "SLWgt0", "SLWeqNA")

#-------------------------------------------------------
# dmax vs FZDZ max ints: derived statistics
#-------------------------------------------------------
#   for FZDZ
#FRAC.liq.FZDZ            <- round((num.dmax.gt00.FZDZ.gt.00 + num.dmax.gt00.FZDZ.eq.NA) / num.total, digits=2)
#FRAC.liq                 <- round(num.dmax.gt00 / num.total, digits=2)
FRAC.liq                 <- round(num.dmax.gt00.cdp.gt01 / num.total, digits=2)
FRAC.AppO                <- round(num.ac.appO / num.total, digits=2)
#POD.Y.hit.FZDZ           <- round(num.dmax.gt00.FZDZ.gt.60 / (num.dmax.gt00.FZDZ.gt.00 + num.dmax.gt00.FZDZ.eq.NA), digits=2)
POD.Y.hit.FZDZ           <- round(num.FZDZ.notNA / num.ac.appO, digits=2)
FZDZ.med.ac.appO         <- median(rawData2$FZDZ.int.match.matrix.max[ind.ac.appO], na.rm=TRUE)

#POD.N.FZDZ               <- round(num.dmax.eq00.FZDZ.eq.NA / (num.dmax.eq00.FZDZ.gt.00 + num.dmax.eq00.FZDZ.eq.NA), digits=2)
POD.N.FZDZ               <- round(num.dmax.eq00.FZDZ.eqNA / num.dmax.eq00, digits=2)
FAR.FZDZ                 <- round(num.dmax.eq00.FZDZ.gt.00 / (num.dmax.eq00.FZDZ.gt.00 + num.dmax.eq00.FZDZ.eq.NA), digits=2)
#MISSR.FZDZ               <- round(num.dmax.gt00.FZDZ.eq.NA / (num.dmax.gt00.FZDZ.eq.NA + num.dmax.gt00.FZDZ.gt.00), digits=2) 
MISSR.FZDZ               <- round(num.dmax.gt50.FZDZ.eqNA / num.dmax.gt50, digits=2) # frac of total that FZDZ had as NA 

#   for MPHA
FRAC.liq.MPHA            <- round((num.dmax.gt00.MPHA.gt.00 + num.dmax.gt00.MPHA.eq.NA) / num.total, digits=2)
POD.Y.hit.MPHA           <- round(num.dmax.gt00.MPHA.gt.60 / (num.dmax.gt00.MPHA.gt.00 + num.dmax.gt00.MPHA.eq.NA), digits=2)
POD.Y.tot.MPHA           <- round(num.dmax.gt00.MPHA.gt.00 / (num.dmax.gt00.MPHA.gt.00 + num.dmax.gt00.MPHA.eq.NA), digits=2)
MPHA.med.ac.mpha         <- median(rawData2$MIXPHA.int.match.matrix.max[ind.ac.mixpha], na.rm=TRUE)
POD.N.MPHA               <- round(num.dmax.eq00.MPHA.eq.NA / (num.dmax.eq00.MPHA.gt.00 + num.dmax.eq00.MPHA.eq.NA), digits=2)
FAR.MPHA                 <- round(num.dmax.eq00.MPHA.gt.00 / (num.dmax.eq00.MPHA.gt.00 + num.dmax.eq00.MPHA.eq.NA), digits=2)
MISSR.MPHA               <- round(num.dmax.gt00.MPHA.eq.NA / (num.dmax.gt00.MPHA.eq.NA + num.dmax.gt00.MPHA.gt.00), digits=2)

#   for SLW
FRAC.AppC               <- round(num.ac.appC / num.total, digits=2)
FRAC.liq.SLW            <- round((num.dmax.gt00.SLW.gt.00 + num.dmax.gt00.SLW.eq.NA) / num.total, digits=2)
POD.Y.hit.SLW           <- round(num.SLW.notNA / num.ac.appC, digits=2)
SLW.med.ac.appC         <- median(rawData2$SLW.int.match.matrix.max[ind.ac.appC], na.rm=TRUE)
POD.N.SLW               <- round(num.dmax.eq00.SLW.eq.NA / (num.dmax.eq00.SLW.gt.00 + num.dmax.eq00.SLW.eq.NA), digits=2)
FAR.SLW                 <- round(num.dmax.eq00.SLW.gt.00 / (num.dmax.eq00.SLW.gt.00 + num.dmax.eq00.SLW.eq.NA), digits=2)
MISSR.SLW               <- round(num.dmax.gt00.SLW.eq.NA / (num.dmax.gt00.SLW.eq.NA + num.dmax.gt00.SLW.gt.00), digits=2)

#-------------------------------------------------------

#----------------------------------------------
# create 4th degree polynomial fits for RadIA ints
#----------------------------------------------
# for dmax v FZDZ int
x      <- rawData2$Dmax.um[ind.ac.appO]
y      <- rawData2$FZDZ.int.match.column.max[ind.ac.appO]
if (length(y) > 0) {
  fit    <- lm(y ~ poly(x, 1, raw = TRUE))
}
xx     <- seq(0, max(rawData2$Dmax.um, na.rm=TRUE)*0.75, length=841)

# for dmax v SLW int
x      <- rawData2$Dmax.um[ind.ac.appC]
y2     <- rawData2$SLW.int.match.column.max[ind.ac.appC]
if (length(y2) > 0) {
  fit2   <- lm(y2 ~ poly(x, 1, raw = TRUE))
}
xx2    <- seq(50, max(rawData2$Dmax.um, na.rm=TRUE)*0.75, length=841)

# for dmax v MIXPHA int
x      <- rawData2$Dmax.um[ind.ac.mixpha]
y3     <- rawData2$MIXPHA.int.match.matrix.max[ind.ac.mixpha]
if (length(y3) > 0) {
  fit3   <- lm(y3 ~ poly(x, 4, raw = TRUE))
}
xx3     <- seq(50, max(rawData2$Dmax.um, na.rm=TRUE)*0.75, length=841)

# for massliqfrac vs FZDZ int
x        <- rawData2$mass.liq.frac[ind.dmax.gt00.FZDZ.gt.00]
y9       <- rawData2$FZDZ.int.match.matrix.max[ind.dmax.gt00.FZDZ.gt.00]
if (length(y9) > 0) {
  fit9    <- lm(y9 ~ poly(x, 4, raw = TRUE))
}
xx9     <- seq(0, max(rawData2$mass.liq.frac, na.rm=TRUE), length=841)

# for massliqfrac vs MIXPHA int
x        <- rawData2$mass.liq.frac[ind.dmax.gt00.FZDZ.gt.00]
y10      <- rawData2$MIXPHA.int.match.matrix.max[ind.dmax.gt00.FZDZ.gt.00]
if (length(y10) > 0) {
  fit10    <- lm(y10 ~ poly(x, 4, raw = TRUE))
}
xx10     <- seq(0, max(rawData2$mass.liq.frac, na.rm=TRUE), length=841)

# for dist from radar vs fracFZDZisnotNA
range.gates.corr   <- 10
mean.range.gates   <- c( 325, 320, 384, 406, 292, 476, 394, 319)
frac.FZDZ.isnot.NA <- c(21/28, 39/44, 13/62, 8/31, 3/8, 0/12, 0/4, 2/2)
x        <- (mean.range.gates+range.gates.corr)*range.gate.dist.km
y20      <- frac.FZDZ.isnot.NA
if (length(y20) > 0) {
  fit20  <- lm(y20 ~ poly(x, 2, raw=TRUE))
}
xx20     <- seq(0, max(mean.range.gates*range.gate.dist.km, na.rm=TRUE, length=length(mean.range.gates)))

#----------------------------------------------

# set up color palatte
n_palette <- 11
marker    <- list(color = rev(brewer.pal(n_palette, "Spectral")))

#----------------------------------------------
# plot
#----------------------------------------------

# plot ROC curves with AUC values
par(mfrow=c(1, 3))
FAR.FZDZ.order  <- order(FAR.FZDZ.extend)
AUC.FZDZ        <- sum(diff(FAR.FZDZ.extend[FAR.FZDZ.order])*rollmean(PODY.FZDZ.extend[FAR.FZDZ.order], 2))
plot(FAR.FZDZ.extend, PODY.FZDZ.extend, xlim=c(0,1), ylim=c(0,1), pch=21, type="b", lwd=1, xlab="FAR", ylab="PODY", main=paste("FZDZ: AUC = ", round(AUC.FZDZ, digits=2), sep=""))
abline(0, 1, col="red")
points(FAR.FZDZ.extend[16], PODY.FZDZ.extend[16], col="blue", type="p", pch=20, lwd=2)
points(0,1, col="blue", type="p", pch=20, lwd=2)
segments(0,1,FAR.FZDZ.extend[16], PODY.FZDZ.extend[16], lty=2, col="blue")
text(FAR.FZDZ.extend[16]+0.14, PODY.FZDZ.extend[16]-0.05, paste("THR=", thresh.matrix[16], sep=""), col="blue")
text(FAR.FZDZ.extend[16]+0.14, PODY.FZDZ.extend[16]-0.10, paste("N=", length(ind.ac.appO), sep=""), col="blue")
grid()
FAR.SLW.order   <- order(FAR.SLW.extend)
AUC.SLW         <- sum(diff(FAR.SLW.extend[FAR.SLW.order])*rollmean(PODY.SLW.extend[FAR.SLW.order], 2))
plot(FAR.SLW.extend, PODY.SLW.extend, xlim=c(0,1), ylim=c(0,1), pch=21, type="b", lwd=1, xlab="FAR", ylab="PODY", main=paste("SLW: AUC = ", round(AUC.SLW, digits=2), sep=""))
abline(0, 1, col="red")
points(FAR.SLW.extend[40], PODY.SLW.extend[40], col="green", type="p", pch=20, lwd=2)
text(FAR.SLW.extend[40]+0.16, PODY.SLW.extend[40]-0.03, paste("THR=", thresh.matrix[39], sep=""), col="green")
text(FAR.SLW.extend[40]+0.16, PODY.SLW.extend[40]-0.08, paste("N=", length(ind.ac.appC), sep=""), col="green")
points(0,1, col="green", type="p", pch=20, lwd=2)
segments(0,1,FAR.SLW.extend[40], PODY.SLW.extend[40], lty=2, col="green")
grid()
# for 0.2 < MLF < 0.8 & LWC.CDP > 0.02
FAR.MPHA.order   <- order(FAR.MPHA)
AUC.MPHA         <- sum(diff(FAR.MPHA[FAR.MPHA.order])*rollmean(PODY.MPHA[FAR.MPHA.order], 2))
plot(FAR.MPHA, PODY.MPHA, xlim=c(0,1), ylim=c(0,1), pch=21, type="b", lwd=1, xlab="FAR", ylab="PODY", main=paste("MIXPHA: AUC = ", round(AUC.MPHA, digits=2), ", LWC.CDP>=0.05", sep=""))
abline(0, 1, col="red")
points(FAR.MPHA[22], PODY.MPHA[22], col="orange", type="p", pch=20, lwd=2)
text(FAR.MPHA[22]+0.10, PODY.MPHA[22]-0.05, paste("THR=", thresh.matrix[22], sep=""), col="orange")
text(FAR.MPHA[22]+0.10, PODY.MPHA[22]-0.10, paste("N=", length(ind.ac.mixpha), sep=""), col="orange")
points(0,1, col="orange", type="p", pch=20, lwd=2)
segments(0,1,FAR.MPHA[22], PODY.MPHA[22], lty=2, col="orange")
grid()

# plot 0A:
# Question: how far away from radar does RadIA work for?
#   Answer: below...
range.gate.corr      <- 10
dates                <- c( "1/22",  "1/31",   "1/8", "1/18", "2/21", "1/19", "1/21", "1/9")
mean.range.gates     <- c(   325,     320,     384,    406,    292,    476,    394,  319)
dist.to.radar        <- (mean.range.gates+range.gate.corr)*range.gate.dist.km
frac.FZDZ.isnot.NA   <- c( 21/28,   39/44,   13/62,   8/31,    3/8,   0/12,    0/4,  2/2)
tot.30s.periods      <- c(    28,      44,      62,     31,      8,     12,      4,    3)
howfarradar.df       <- data.frame(dates, dist.to.radar, frac.FZDZ.isnot.NA, tot.30s.periods)

plot(frac.FZDZ.isnot.NA~dist.to.radar, data=howfarradar.df, xlim=c(50,140), ylim=c(0,1), pch=20, xlab="range from radar [km]", ylab="RadIA fraction not 'NA'", main="")
with(howfarradar.df, text(frac.FZDZ.isnot.NA~dist.to.radar, labels=paste(dates, ", N=", tot.30s.periods, sep=""), pos=4))
#lines(xx20, predict(fit20, data.frame(x=xx20)), col="grey", lwd=2)
grid()

#-------------------------------------
#plot #1a: dmax from liq.gt10um versus FZDZ 
plot.new()
par(mfrow=c(1,1))
par(mar=c(7, 4, 2, 1.5))
#colPal        <- colorRampPalette(palette(marker$color))
#rawData2$Col4 <- colPal(11)[as.numeric(cut(rawData2$MIXPHA.int.match.matrix.max, breaks=110))]

plot(rawData2$Dmax.um[ind.ac.appO], rawData2$FZDZ.int.match.matrix.max[ind.ac.appO], pch=20, col="blue", lwd=1.0, xlim=c(0, 1000), ylim=c(0,1), xaxt="n", xlab="D-max [um]", ylab="FZDZ int")
#plot(rawData2$Dmax.um[ind.ac.appO], rawData2$FZDZ.int.match.matrix.max[ind.ac.appO], pch=20, col="blue", lwd=1.0, xlim=c(0, 1000), ylim=c(0,1), main="Dmax>100, LWCRos-Nev>0.05, MLF>0.8", xaxt="n", xlab="D-max [um]", ylab="FZDZ int")
# for sized points, use cex=exp(rawData2$LAC.ros[ind.ac.appO])^2.5 in plot command above, with pch=21
#if (length(y) > 0) {
#  lines(xx, predict(fit, data.frame(x=xx)), col="blue", lwd=2)
#}
#legend(250, 0.10, legend=c("FZDZ", "LWC-Ros<0.10", "LWC-Ros=0.35", "LWC-Ros=0.70"), pch=c(20, 21, 21, 21), pt.cex=c(0.001, 1.0, 2.0, 4.0), col=c("blue", "blue", "blue", "blue"), lty=c(1,0,0,0), cex=0.8)
x.ticks <- c(0, 5, 50, 100, 150, 200, 300, 400, 500, 1000, 1500)
axis(side=1, at=x.ticks)
axis(side=2)
abline(h=seq(0, 1, 0.2), v=x.ticks, col="lightgray", lty=3)
abline(h=0.30, lty=2)
#abline(h=0.6, lty=2)
box()
text(50, 0.32, "th=0.30")
#text(50, 0.62, "th=0.60")
text(350, 0.25, paste("FRAC-liq= ",   FRAC.liq,         sep=""), col="blue")
text(350, 0.20, paste("FRAC-AppO= ",  FRAC.AppO,        sep=""), col="blue")
#text(350, 0.17, paste("N-AppO= ",     num.ac.appO,      sep=""), col="blue")
#text(350, 0.14, paste("POD-Y-AppO= ", POD.Y.hit.FZDZ,   sep=""), col="blue")
#text(350, 0.11, paste("FZDZ-med= ",   FZDZ.med.ac.appO, sep=""), col="blue")
#text(350, 0.08, paste("MISSR= ",      MISSR.FZDZ,       sep=""), col="blue")
#text(350, 0.05, paste("POD-N= ",      POD.N.FZDZ,       sep=""), col="blue")
#text(350, 0.02, paste("FAR= ",        FAR.FZDZ,         sep=""), col="blue")

#LWC.nev.bias        <- 0.03
#LWCdiff             <- rawData2$LWC.ros[ind.ac.appO] - rawData2$LWC.nev[ind.ac.appO] + LWC.nev.bias
#ind.LWCdiff.gtneg05 <- which(LWCdiff > -0.05)
#plot(LWCdiff[ind.LWCdiff.gtneg05], rawData2$FZDZ.int.match.matrix.max[ind.LWCdiff.gtneg05], pch=20, col="blue", xlim=c(-0.2, 0.5), ylim=c(0,1), xlab="LWC-Rose/Nev [gm-3]", ylab="RadIA int")
#abline(h=0.6, lty=2)
#abline(v=0.04, lty=2)
#abline(v=-0.04, lty=2)
#text(0.30, 0.10, paste(as.character(num.total), " 30-s averages", sep=""))
#grid()

# plot #1b: SLW when AC-indicated App O
plot(rawData2$SLW.int.match.matrix.max[ind.ac.appO], rawData2$FZDZ.int.match.matrix.max[ind.ac.appO], pch=20, col="blue", lwd=1.0, xlim=c(0, 1), ylim=c(0,1), main="Dmax>100, LWCRos-Nev>0.05, MLF>0.8", xlab="SLW int", ylab="FZDZ int")
abline(v=0.6, lty=2)
abline(h=0.6, lty=2)
text(0.20, 0.20, paste("FZDZint-med@AppO-AC = ", round(median(rawData2$FZDZ.int.match.matrix.max[ind.ac.appO], na.rm=TRUE), digits=2), sep=""), col="blue")
text(0.20, 0.15, paste("SLWint-med@AppO-AC = ",  round(median(rawData2$SLW.int.match.matrix.max[ind.ac.appO],  na.rm=TRUE), digits=2), sep=""), col="blue")
text(0.20, 0.10, paste("N = ", length(!is.na(rawData2$FZDZ.int.match.matrix.max[ind.ac.appO])), sep=""), col="blue")
grid()

# plot #1c: MIXPHA when AC-indicated App O
plot(rawData2$MIXPHA.int.match.matrix.max[ind.ac.appO], rawData2$FZDZ.int.match.matrix.max[ind.ac.appO], pch=20, col="blue", lwd=1.0, xlim=c(0, 1), ylim=c(0,1), main="Dmax>100, LWCRos-Nev>0.05, MLF>0.8", xlab="MIXPHA int", ylab="FZDZ int")
abline(v=0.6, lty=2)
abline(h=0.6, lty=2)
text(0.20, 0.20, paste("FZDZint-med@AppO-AC = ", round(median(rawData2$FZDZ.int.match.matrix.max[ind.ac.appO],   na.rm=TRUE), digits=2), sep=""), col="blue")
text(0.20, 0.15, paste("MPHAint-med@AppO-AC = ", round(median(rawData2$MIXPHA.int.match.matrix.max[ind.ac.appO], na.rm=TRUE), digits=2), sep=""), col="blue")
text(0.20, 0.10, paste("N = ", length(ind.ac.appO), sep=""), col="blue")
grid()

#-----------
#plot #2a: dmax from liq.gt10um versus MIXPHA
plot.new()
par(mfrow=c(1,1))
par(mar=c(7, 4, 2, 1.5))
#colPal        <- colorRampPalette(palette(marker$color))
#rawData2$Col3 <- colPal(11)[as.numeric(cut(rawData2$mass.liq.frac[ind.ac.mixpha], breaks = 9))]
plot(rawData2$Dmax.um[ind.ac.mixpha], rawData2$MIXPHA.int.match.matrix.max[ind.ac.mixpha], pch=20, col="orange", lwd=1.0, xlim=c(0, 1000), ylim=c(0,1), xaxt="n", xlab="D-max [um]", ylab="MIXPHA int")
#plot(rawData2$Dmax.um[ind.ac.mixpha], rawData2$MIXPHA.int.match.matrix.max[ind.ac.mixpha], pch=20, col="orange", lwd=1.0, xlim=c(0, 1000), ylim=c(0,1), main="CDP>0.01, 0.20<MLF<0.80, column", xaxt="n", xlab="D-max [um]", ylab="MIXPHA int")
# cex=exp(rawData2$mass.liq.frac[ind.ac.appC])^1.2 with pch=21
#set.colorbar(cbpos="b", pal=marker$color, zlim=c(0.00, max(rawData2$mass.liq.frac[ind.dmax.gt00], na.rm=TRUE)))
lines(xx3, predict(fit3, data.frame(x=xx3)), col="orange", lwd=2)
#legend(250, 0.10, legend=c("MPHA", "MLF~0.0", "MLF=0.5", "MLF~1.0"), pch=c(20, 21, 21, 21), pt.cex=c(0.01, 1.0, 2.0, 3.0), col=c("orange", "orange", "orange", "orange"), lty=c(1,0,0,0), cex=0.8)
x.ticks <- c(0, 5, 50, 100, 150, 200, 300, 400, 500, 1000, 1500)
axis(side=1, at=x.ticks)
axis(side=2)
abline(h=seq(0, 1, 0.2), v=x.ticks, col="lightgray", lty=3)
abline(h=0.42, lty=2)
box()
text(50, 0.44, "th=0.42")
text(400, 0.25, paste("FRAC-liq= ",  FRAC.liq,       sep=""),             col="orange")
text(400, 0.20, paste("FRAC-MPHA= ", round(293/num.total, digits=2), sep=""), col="orange")
#text(400, 0.17, "N-MIXPHA= 293",                                          col="orange")
#text(400, 0.14, paste("POD-Y-MIXPHA= ", round(168/293, digits=2), sep=""),                 col="orange")
#text(400, 0.11, paste("MIXPHA-med= ", MPHA.med.ac.mpha, sep=""), col="orange")
#text(400, 0.10, paste("MISSR= ",     MISSR.MPHA,     sep=""), col="orange")
#text(400, 0.07, paste("POD-N= ",     POD.N.MPHA,     sep=""), col="orange")
#text(400, 0.04, paste("FAR= ",       FAR.MPHA,       sep=""), col="orange")

LWC.ros.bias <- 0.03
plot(rawData2$LWC.ros[ind.dmax.gt00.FZDZ.gt.00] - rawData2$LWC.nev[ind.dmax.gt00.FZDZ.gt.00] + LWC.ros.bias, rawData2$MIXPHA.int.match.matrix.max[ind.dmax.gt00], pch=20, col="orange", xlim=c(-0.2, 0.5), ylim=c(0,1), xlab="LWC-Nev/CDP [gm-3]", ylab="RadIA int")
abline(h=0.6, lty=2)
abline(v=0.04, lty=2)
abline(v=-0.04, lty=2)
text(0.30, 0.10, paste(as.character(num.total), " 30-s averages", sep=""))
grid()

# plot #2b: FZDZ when AC-indicated MIXPHA
plot(rawData2$FZDZ.int.match.matrix.max[ind.ac.mixpha], rawData2$MIXPHA.int.match.matrix.max[ind.ac.mixpha], pch=20, col="orange", lwd=1.0, xlim=c(0, 1), ylim=c(0,1), main="CDP>0.01, 0.20<MLF<0.80", xlab="FZDZ int", ylab="MIXPHA int")
abline(v=0.6, lty=2)
abline(h=0.6, lty=2)
text(0.20, 0.20, paste("MIXPHAint-med@MIXPHA-AC = ", round(median(rawData2$MIXPHA.int.match.matrix.max[ind.ac.mixpha], na.rm=TRUE), digits=2), sep=""), col="orange")
text(0.20, 0.15, paste("FZDZint-med@MIXPHA-AC = ",    round(median(rawData2$FZDZ.int.match.matrix.max[ind.ac.mixpha],  na.rm=TRUE), digits=2), sep=""), col="orange")
text(0.20, 0.10, paste("N = ", length(ind.ac.mixpha), sep=""), col="orange")
grid()

# plot #2c: SLW when AC-indicated MIXPHA
plot(rawData2$SLW.int.match.matrix.max[ind.ac.mixpha], rawData2$MIXPHA.int.match.matrix.max[ind.ac.mixpha], pch=20, col="orange", lwd=1.0, xlim=c(0, 1), ylim=c(0,1), main="CDP>0.01, 0.20<MLF<0.80", xlab="SLW int", ylab="MIXPHA int")
abline(v=0.6, lty=2)
abline(h=0.6, lty=2)
text(0.20, 0.20, paste("MIXPHAint-med@MIXPHA-AC = ", round(median(rawData2$MIXPHA.int.match.matrix.max[ind.ac.mixpha],   na.rm=TRUE), digits=2), sep=""), col="orange")
text(0.20, 0.15, paste("SLWint-med@MIXPHA-AC = ",    round(median(rawData2$SLW.int.match.matrix.max[ind.ac.mixpha],      na.rm=TRUE), digits=2), sep=""), col="orange")
text(0.20, 0.10, paste("N = ", length(ind.ac.mixpha), sep=""), col="orange")
grid()

#-----------
#plot #3: dmax from liq.gt10um versus SLW
plot.new()
par(mfrow=c(1,1))
par(mar=c(7, 4, 2, 1.5))
#colPal        <- colorRampPalette(palette(marker$color))
#rawData2$Col5 <- colPal(11)[as.numeric(cut(rawData2$LWC.nev[ind.dmax.gt00], breaks=110))]
plot(rawData2$Dmax.um[ind.ac.appC], rawData2$SLW.int.match.column.max[ind.ac.appC], pch=20, col="green", lwd=1.0, xlim=c(0, 1000), ylim=c(0,1), xaxt="n", xlab="D-max [um]", ylab="SLW int")
#plot(rawData2$Dmax.um[ind.ac.appC], rawData2$SLW.int.match.column.max[ind.ac.appC], pch=20, col="green", lwd=1.0, xlim=c(0, 1000), ylim=c(0,1), main="Dmax<100, LWCRos-Nev<0.05, MLF>0.2, column", xaxt="n", xlab="D-max [um]", ylab="SLW int")
#cex=exp(3*rawData2$LWC.nev[ind.ac.appC]), with pch=21
lines(xx2, predict(fit2, data.frame(x=xx2)), col="green", lwd=2)
#legend(250, 0.10, legend=c("SLW", "LWC-NEV<0.10", "LWC-NEV=0.35", "LWC-NEV=0.60"), pch=c(20, 21, 21, 21), pt.cex=c(0.01, 1.3, 2.0, 3.0), col=c("green", "green", "green", "green"), lty=c(1,0,0,0), cex=0.8)
x.ticks <- c(0, 5, 50, 100, 150, 200, 300, 400, 500, 1000, 1500)
axis(side=1, at=x.ticks)
axis(side=2)
abline(h=seq(0, 1, 0.2), v=x.ticks, col="lightgray", lty=3)
abline(h=0.76, lty=2)
box()
text(400, 0.25, paste("FRAC-liq= ",   FRAC.liq,        sep=""), col="green")
text(400, 0.20, paste("FRAC-AppC= ",  FRAC.AppC,       sep=""), col="green")
#text(400, 0.17, paste("N-AppC= ",     num.ac.appC,     sep=""), col="green")
#text(400, 0.14, paste("POD-Y-AppC= ", POD.Y.hit.SLW,   sep=""), col="green")
#text(400, 0.11, paste("SLW-med= ",    SLW.med.ac.appC, sep=""), col="green")
#text(400, 0.11, paste("MISSR= ",     MISSR.SLW,        sep=""), col="green")
#text(400, 0.08, paste("POD-N= ",     POD.N.SLW,        sep=""), col="green")
#text(400, 0.05, paste("FAR= ",       FAR.SLW,          sep=""), col="green")
text(250, 0.78, "th=0.76")

#LWC.ros.bias <- 0.03
#plot(rawData2$LWC.ros[ind.dmax.gt00] - rawData2$LWC.nev[ind.dmax.gt00] + LWC.ros.bias, rawData2$SLW.int.match.matrix.max[ind.dmax.gt00], pch=20, col="blue", xlim=c(-0.2, 0.5), ylim=c(0,1), xlab="LWC-Nev/CDP [gm-3]", ylab="RadIA int")
#abline(h=0.6, lty=2)
#abline(v=0.04, lty=2)
#abline(v=-0.04, lty=2)
#text(0.30, 0.10, paste(as.character(num.total), " 30-s averages", sep=""))
#grid()

# plot #3b: FZDZ int when AC-indicated App C
plot(rawData2$FZDZ.int.match.matrix.max[ind.ac.appC], rawData2$SLW.int.match.matrix.max[ind.ac.appC], pch=20, col="green", lwd=1.0, xlim=c(0, 1), ylim=c(0,1), main="Dmax<100,LWCRos-Nev<0.05,MLF>0.2", xlab="FZDZ int", ylab="SLW int")
abline(v=0.6, lty=2)
abline(h=0.6, lty=2)
text(0.20, 0.20, paste("SLWint-med@AppC-AC = ",   round(median(rawData2$SLW.int.match.matrix.max[ind.ac.appC],   na.rm=TRUE), digits=2), sep=""), col="green")
text(0.20, 0.15, paste("FZDZint-med@AppC-AC = ",  round(median(rawData2$FZDZ.int.match.matrix.max[ind.ac.appC],  na.rm=TRUE), digits=2), sep=""), col="green")
text(0.20, 0.10, paste("N = ", length(ind.ac.appC), sep=""), col="green")
grid()

# plot #3c: MIXPHA int when AC-indicated APP C
plot(rawData2$MIXPHA.int.match.matrix.max[ind.ac.appC], rawData2$SLW.int.match.matrix.max[ind.ac.appC], pch=20, col="green", lwd=1.0, xlim=c(0, 1), ylim=c(0,1), main="Dmax<100,LWCRos-Nev<0.05,MLF>0.2", xlab="MIXPHA int", ylab="SLW int")
abline(v=0.6, lty=2)
abline(h=0.6, lty=2)
text(0.20, 0.20, paste("SLWint-med@AppC-AC = ",    round(median(rawData2$SLW.int.match.matrix.max[ind.ac.appC],   na.rm=TRUE), digits=2), sep=""), col="green")
text(0.20, 0.15, paste("MIXPHAint-med@AppC-AC = ", round(median(rawData2$MIXPHA.int.match.matrix.max[ind.ac.appC],      na.rm=TRUE), digits=2), sep=""), col="green")
text(0.20, 0.10, paste("N = ", length(ind.ac.appC), sep=""), col="green")
grid()

#-------------------------
#plot #XXX: dmax for ICE ONLY (3xLWC == 0 & mass.liq.frac <= 0.05)
plot.new()
par(mfrow=c(2,2))
par(mar=c(7, 4, 2, 1.5))
#colPal              <- colorRampPalette(palette(marker$color))
#rawData2$Col5       <- colPal(11)[as.numeric(cut(rawData2$LWC.nev[ind.dmax.gt00], breaks=110))]
plot(rawData2$Dmax.um[ind.ac.iceonly3], rawData2$MIXPHA.int.match.column.max[ind.ac.iceonly3], pch=20, col="orange", lwd=1.0, xlim=c(0, 1000), ylim=c(0,1), main="LWC=0, Dmax-liq=0", xaxt="n", xlab="D-max [um]", ylab="MIXPHA int")
#cex=exp(3*rawData2$LWC.nev[ind.ac.appC]), with pch=21
#lines(xx2, predict(fit2, data.frame(x=xx2)), col="green", lwd=2)
#legend(250, 0.10, legend=c("SLW", "LWC-NEV<0.10", "LWC-NEV=0.35", "LWC-NEV=0.60"), pch=c(20, 21, 21, 21), pt.cex=c(0.01, 1.3, 2.0, 3.0), col=c("green", "green", "green", "green"), lty=c(1,0,0,0), cex=0.8)
x.ticks <- c(0, 5, 50, 100, 150, 200, 300, 400, 500, 1000, 1500)
axis(side=1, at=x.ticks)
axis(side=2)
abline(h=seq(0, 1, 0.2), v=x.ticks, col="lightgray", lty=3)
abline(h=0.60, lty=2)
text(50, 0.62, "th=0.60")
text(400, 0.15, paste("MIXPHA-mean= ",   round(mean(rawData2$MIXPHA.int.match.column.max[ind.ac.iceonly3]), digits=2),  sep=""), col="orange")
text(400, 0.11, paste("MIXPHA-med= ",    median(rawData2$MIXPHA.int.match.column.max[ind.ac.iceonly3]),                 sep=""), col="orange")
text(400, 0.07, paste("N= ",             length(ind.ac.iceonly3),                                                       sep=""), col="orange")
plot(rawData2$Dmax.um[ind.ac.iceonly3], rawData2$FZDZ.int.match.column.max[ind.ac.iceonly3], pch=20, col="blue", lwd=1.0, xlim=c(0, 1000), ylim=c(0,1), main="LWC=0, Dmax-liq=0", xaxt="n", xlab="D-max [um]", ylab="FZDZ int")
#cex=exp(3*rawData2$LWC.nev[ind.ac.appC]), with pch=21
#lines(xx2, predict(fit2, data.frame(x=xx2)), col="green", lwd=2)
#legend(250, 0.10, legend=c("SLW", "LWC-NEV<0.10", "LWC-NEV=0.35", "LWC-NEV=0.60"), pch=c(20, 21, 21, 21), pt.cex=c(0.01, 1.3, 2.0, 3.0), col=c("green", "green", "green", "green"), lty=c(1,0,0,0), cex=0.8)
x.ticks <- c(0, 5, 50, 100, 150, 200, 300, 400, 500, 1000, 1500)
axis(side=1, at=x.ticks)
axis(side=2)
abline(h=seq(0, 1, 0.2), v=x.ticks, col="lightgray", lty=3)
abline(h=0.60, lty=2)
text(100, 0.62, "th=0.60")
text(400, 0.15, paste("FZDZ-mean= ",   round(mean(rawData2$FZDZ.int.match.column.max[ind.ac.iceonly3], na.rm=TRUE), digits=2),  sep=""), col="blue")
text(400, 0.11, paste("FZDZ-med= ",    median(rawData2$FZDZ.int.match.column.max[ind.ac.iceonly3],     na.rm=TRUE),             sep=""), col="blue")
text(400, 0.07, paste("N= ",           length(ind.ac.iceonly3),                                                                 sep=""), col="blue")
plot(rawData2$Dmax.um[ind.ac.iceonly3], rawData2$SLW.int.match.column.max[ind.ac.iceonly3], pch=20, col="green", lwd=1.0, xlim=c(0, 1000), ylim=c(0,1), main="LWC=0, Dmax-liq=0", xaxt="n", xlab="D-max [um]", ylab="SLW int")
#cex=exp(3*rawData2$LWC.nev[ind.ac.appC]), with pch=21
#lines(xx2, predict(fit2, data.frame(x=xx2)), col="green", lwd=2)
#legend(250, 0.10, legend=c("SLW", "LWC-NEV<0.10", "LWC-NEV=0.35", "LWC-NEV=0.60"), pch=c(20, 21, 21, 21), pt.cex=c(0.01, 1.3, 2.0, 3.0), col=c("green", "green", "green", "green"), lty=c(1,0,0,0), cex=0.8)
x.ticks <- c(0, 5, 50, 100, 150, 200, 300, 400, 500, 1000, 1500)
axis(side=1, at=x.ticks)
axis(side=2)
abline(h=seq(0, 1, 0.2), v=x.ticks, col="lightgray", lty=3)
#abline(h=0.60, lty=2)
abline(h=0.75, lty=2)
text(100, 0.77, "th=0.75")
text(400, 0.15, paste("SLW-mean= ",   round(mean(rawData2$SLW.int.match.column.max[ind.ac.iceonly3], na.rm=TRUE), digits=2),  sep=""), col="green")
text(400, 0.11, paste("SLW-med= ",    median(rawData2$SLW.int.match.column.max[ind.ac.iceonly3],     na.rm=TRUE),             sep=""), col="green")
text(400, 0.07, paste("N= ",          length(ind.ac.iceonly3),                                                                sep=""), col="green")

#-----------
# plot #4 : FZDZ vs mass liq frac
plot.new()
par(mfrow=c(1,1))
par(mar=c(11,4,4,4))
colPal       <- colorRampPalette(palette(marker$color))
rawData2$Col <- colPal(11)[as.numeric(cut(rawData2$LWC.ros, breaks = 110))]
plot(rawData2$mass.liq.frac[ind.dmax.gt00.FZDZ.gt.00], rawData2$FZDZ.int.match.matrix.max[ind.dmax.gt00.FZDZ.gt.00], pch=20, col=rawData2$Col, xlim=c(0, 1), ylim=c(0,1), xlab="mass liquid fraction", ylab="RadIA FZDZ int")
#plot(rawData2$mass.liq.frac[ind.dmax.gt00.FZDZ.gt.00], rawData2$FZDZ.int.match.matrix.max[ind.dmax.gt00.FZDZ.gt.00], pch=20, col="blue", xlim=c(0, 1), ylim=c(0,1), xlab="mass liquid fraction", ylab="RadIA FZDZ int")
set.colorbar(cbpos="b", pal=marker$color, zlim=c(0.00, max(rawData2$LWC.ros[ind.dmax.gt00.FZDZ.gt.00], na.rm=TRUE)))
if (length(y10) > 0) {
  lines(xx9, predict(fit9, data.frame(x=xx9)), col="blue", lwd=2)
}
abline(h=0.60, lty=2)
abline(h=0.65, lty=2)
abline(v=0.30, lty=2)
abline(v=0.70, lty=2)
#segments( 0.7,  -0.10, 0.7, 0.65, lty=2, lwd=3)
#segments( 0.7,   0.65, 1.1, 0.65, lty=2, lwd=3)
#segments( 0.3,   1.10, 0.3, 0.60, lty=2, lwd=3)
#segments(-0.1,   0.60, 0.3, 0.60, lty=2, lwd=3)
legend(0.5, 0.20, legend=c("FZDZ"), col=c("blue"), lty=1, cex=0.8)
text(0.5, 0.10, paste("Data represents ", as.character(num.total), " 30-s averages ", "(", as.character(round(num.total/2/60), digits=1), " hrs)", sep=""))
text(0.5, 0.07, paste("over ", as.character(num.cases), " case dates", sep=""))
grid()

#ind.mlf.gt.70.MPHA.lt.65 <- which(rawData2$mass.liq.frac[ind.dmax.gt00.FZDZ.gt.00] > 0.70 & rawData2$MIXPHA.int.match.matrix.max[ind.dmax.gt00.FZDZ.gt.00] <= 0.65)
#ind.mlf.gt.70.MPHA.gt.65 <- which(rawData2$mass.liq.frac[ind.dmax.gt00.FZDZ.gt.00] > 0.70 & rawData2$MIXPHA.int.match.matrix.max[ind.dmax.gt00.FZDZ.gt.00] > 0.65)
#ind.mlf.lt.30.MPHA.gt.60 <- which(rawData2$mass.liq.frac[ind.dmax.gt00.FZDZ.gt.00] < 0.30 & rawData2$MIXPHA.int.match.matrix.max[ind.dmax.gt00.FZDZ.gt.00] >= 0.60)
#ind.mlf.lt.30.MPHA.lt.60 <- which(rawData2$mass.liq.frac[ind.dmax.gt00.FZDZ.gt.00] < 0.30 & rawData2$MIXPHA.int.match.matrix.max[ind.dmax.gt00.FZDZ.gt.00] < 0.60)
#frac.mlf.gt.70.MPHA.lt.65 <- length(ind.mlf.gt.70.MPHA.lt.65) / (length(ind.mlf.gt.70.MPHA.lt.65) + length(ind.mlf.gt.70.MPHA.gt.65))
#frac.mlf.lt.30.MPHA.gt.60 <- length(ind.mlf.lt.30.MPHA.gt.60) / (length(ind.mlf.lt.30.MPHA.gt.60) + length(ind.mlf.lt.30.MPHA.lt.60))
#text(0.12, 0.95, paste(round(frac.mlf.lt.30.MPHA.gt.60, digits=2)*100, "% pts MLF<0.3 over 0.60 int", sep=""))
#text(0.88, 0.40, paste(round(frac.mlf.gt.70.MPHA.lt.65, digits=2)*100, "% pts MLF>0.7 under 0.65 int", sep=""))

#-----------
# plot #5 : MIXPHA vs mass liq frac, colored by LWC-ROS
plot.new()
par(mfrow=c(1,1))
par(mar=c(11,4,4,4))
colPal        <- colorRampPalette(palette(marker$color))
rawData2$Col2 <- colPal(11)[as.numeric(cut(rawData2$LWC.ros, breaks = 110))]
#plot(rawData2$mass.liq.frac[ind.dmax.gt00.FZDZ.gt.00], rawData2$MIXPHA.int.match.matrix.max[ind.dmax.gt00.FZDZ.gt.00], pch=20, col="orange", xlim=c(0, 1), ylim=c(0,1), xlab="mass liquid fraction", ylab="RadIA MIXPHA int")
plot(rawData2$mass.liq.frac[ind.dmax.gt00.FZDZ.gt.00], rawData2$MIXPHA.int.match.matrix.max[ind.dmax.gt00.FZDZ.gt.00], pch=20, col=rawData2$Col2[ind.dmax.gt00.FZDZ.gt.00], xlim=c(0, 1), ylim=c(0,1), xlab="mass liquid fraction", ylab="RadIA MIXPHA int")
#set.colorbar(cbpos="b", pal=marker$color, zlim=c(0.00, max(rawData2$LWC.nev[ind.dmax.gt00.FZDZ.gt.00], na.rm=TRUE)))
set.colorbar(cbpos="b", pal=marker$color, zlim=c(0.00, max(rawData2$LWC.ros[ind.dmax.gt00.FZDZ.gt.00], na.rm=TRUE)))
if (length(y10) > 0) {
  lines(xx10, predict(fit10, data.frame(x=xx10)), col="orange", lwd=2)
}
abline(h=0.60, lty=2)
abline(h=0.65, lty=2)
abline(v=0.33, lty=2)
abline(v=0.66, lty=2)
segments(0.66, -0.1, 0.66, 0.65, lty=2, lwd=3)
segments(0.66, 0.65, 1.10, 0.65, lty=2, lwd=3)
segments(0.33, 1.10, 0.33, 0.60, lty=2, lwd=3)
segments(-0.1, 0.60, 0.33, 0.60, lty=2, lwd=3)
legend(0.45, 0.30, legend=c("MPHA"), col=c("orange"), lty=1, cex=0.8)
text(  0.50, 0.20, paste(as.character(length(ind.dmax.gt00.FZDZ.gt.00)), " / ", as.character(num.total.pts), " 30-s averages ", "(", as.character(round(num.total/2/60), digits=1), " hrs)", sep=""))
text(  0.50, 0.15, paste("over ", as.character(num.cases), " case dates", sep=""))
grid()

ind.mlf.lt.00               <- which(rawData2$mass.liq.frac[ind.dmax.gt00.FZDZ.gt.00] < 0.00)
ind.mlf.gt.66.MPHA.lt.65    <- which(rawData2$mass.liq.frac[ind.dmax.gt00.FZDZ.gt.00] > 0.66 & rawData2$MIXPHA.int.match.matrix.max[ind.dmax.gt00.FZDZ.gt.00] <= 0.65)
ind.mlf.gt.66.MPHA.gt.65    <- which(rawData2$mass.liq.frac[ind.dmax.gt00.FZDZ.gt.00] > 0.66 & rawData2$MIXPHA.int.match.matrix.max[ind.dmax.gt00.FZDZ.gt.00] > 0.65)
ind.mlf.lt.33.MPHA.gt.60    <- which(rawData2$mass.liq.frac[ind.dmax.gt00.FZDZ.gt.00] < 0.33 & rawData2$MIXPHA.int.match.matrix.max[ind.dmax.gt00.FZDZ.gt.00] >= 0.60)
ind.mlf.lt.33.MPHA.lt.60    <- which(rawData2$mass.liq.frac[ind.dmax.gt00.FZDZ.gt.00] < 0.33 & rawData2$MIXPHA.int.match.matrix.max[ind.dmax.gt00.FZDZ.gt.00] < 0.60)
frac.mlf.gt.66.MPHA.lt.65   <- length(ind.mlf.gt.66.MPHA.lt.65) / (length(ind.mlf.gt.66.MPHA.lt.65) + length(ind.mlf.gt.66.MPHA.gt.65))
frac.mlf.lt.33.MPHA.gt.60   <- length(ind.mlf.lt.33.MPHA.gt.60) / (length(ind.mlf.lt.33.MPHA.gt.60) + length(ind.mlf.lt.33.MPHA.lt.60))
text(0.15, 1.02, paste(round(frac.mlf.lt.33.MPHA.gt.60, digits=2)*100, "% pts MLF<0.33 over 0.60 int", sep=""))
text(0.85, 0.30, paste(round(frac.mlf.gt.66.MPHA.lt.65, digits=2)*100, "% pts MLF>0.66 under 0.65 int", sep=""))

ind.mlf.gt.66               <- which(rawData2$mass.liq.frac[ind.dmax.gt00.FZDZ.gt.00] > 0.66)
ind.mlf.bt.33.66            <- which(rawData2$mass.liq.frac[ind.dmax.gt00.FZDZ.gt.00] >= 0.33 & rawData2$mass.liq.frac[ind.dmax.gt00.FZDZ.gt.00] <= 0.66)
ind.mlf.lt.33               <- which(rawData2$mass.liq.frac[ind.dmax.gt00.FZDZ.gt.00] < 0.33)
LWC.ros.median.mlf.gt.66    <- round(median(rawData2$LWC.ros[ind.dmax.gt00.FZDZ.gt.00[ind.mlf.gt.66]], na.rm=TRUE),    digits=2)
LWC.ros.median.mlf.bt.33.66 <- round(median(rawData2$LWC.ros[ind.dmax.gt00.FZDZ.gt.00[ind.mlf.bt.33.66]], na.rm=TRUE), digits=2)
LWC.ros.median.mlf.lt.33    <- round(median(rawData2$LWC.ros[ind.dmax.gt00.FZDZ.gt.00[ind.mlf.lt.33]], na.rm=TRUE),    digits=2)
LWC.ros.mean.mlf.gt.66      <- round(mean(rawData2$LWC.ros[ind.dmax.gt00.FZDZ.gt.00[ind.mlf.gt.66]], na.rm=TRUE),      digits=2)
LWC.ros.mean.mlf.bt.33.66   <- round(mean(rawData2$LWC.ros[ind.dmax.gt00.FZDZ.gt.00[ind.mlf.bt.33.66]], na.rm=TRUE),   digits=2)
LWC.ros.mean.mlf.lt.33      <- round(mean(rawData2$LWC.ros[ind.dmax.gt00.FZDZ.gt.00[ind.mlf.lt.33]], na.rm=TRUE),      digits=2)
text(0.16, 0.05, paste("LWC-median=", LWC.ros.median.mlf.lt.33,    sep=""))
text(0.50, 0.05, paste("LWC-median=", LWC.ros.median.mlf.bt.33.66, sep=""))
text(0.84, 0.05, paste("LWC-median=", LWC.ros.median.mlf.gt.66,    sep=""))

#---------------
# plot #5b : MIXPHA vs mass liq frac, colored by D-99
plot.new()
par(mfrow=c(1,1))
par(mar=c(11,4,4,4))
colPal        <- colorRampPalette(palette(marker$color))
rawData2$Col2 <- colPal(11)[as.numeric(cut(rawData2$dmax.um, breaks = 110))]
#plot(rawData2$mass.liq.frac[ind.dmax.gt00.FZDZ.gt.00], rawData2$MIXPHA.int.match.matrix.max[ind.dmax.gt00.FZDZ.gt.00], pch=20, col="orange", xlim=c(0, 1), ylim=c(0,1), xlab="mass liquid fraction", ylab="RadIA MIXPHA int")
plot(rawData2$mass.liq.frac[ind.dmax.gt00.FZDZ.gt.00], rawData2$MIXPHA.int.match.matrix.max[ind.dmax.gt00.FZDZ.gt.00], pch=20, col=rawData2$Col2[ind.dmax.gt00.FZDZ.gt.00], xlim=c(0, 1), ylim=c(0,1), xlab="mass liquid fraction", ylab="RadIA MIXPHA int")
set.colorbar(cbpos="b", pal=marker$color, zlim=c(0.00, max(rawData2$dmax.um[ind.dmax.gt00.FZDZ.gt.00], na.rm=TRUE)))
if (length(y10) > 0) {
  lines(xx10, predict(fit10, data.frame(x=xx10)), col="orange", lwd=2)
}
abline(h=0.60, lty=2)
abline(h=0.65, lty=2)
abline(v=0.33, lty=2)
abline(v=0.66, lty=2)
segments(0.66, -0.1, 0.66, 0.65, lty=2, lwd=3)
segments(0.66, 0.65, 1.10, 0.65, lty=2, lwd=3)
segments(0.33, 1.10, 0.33, 0.60, lty=2, lwd=3)
segments(-0.1, 0.60, 0.33, 0.60, lty=2, lwd=3)
legend(0.45, 0.30, legend=c("MPHA"), col=c("orange"), lty=1, cex=0.8)
text(  0.50, 0.20, paste(as.character(length(ind.dmax.gt00.FZDZ.gt.00)), " / ", as.character(num.total.pts), " 30-s averages ", "(", as.character(round(num.total/2/60), digits=1), " hrs)", sep=""))
text(  0.50, 0.15, paste("over ", as.character(num.cases), " case dates", sep=""))
grid()

ind.mlf.lt.00               <- which(rawData2$mass.liq.frac[ind.dmax.gt00.FZDZ.gt.00] < 0.00)
ind.mlf.gt.66.MPHA.lt.65    <- which(rawData2$mass.liq.frac[ind.dmax.gt00.FZDZ.gt.00] > 0.66 & rawData2$MIXPHA.int.match.matrix.max[ind.dmax.gt00.FZDZ.gt.00] <= 0.65)
ind.mlf.gt.66.MPHA.gt.65    <- which(rawData2$mass.liq.frac[ind.dmax.gt00.FZDZ.gt.00] > 0.66 & rawData2$MIXPHA.int.match.matrix.max[ind.dmax.gt00.FZDZ.gt.00] > 0.65)
ind.mlf.lt.33.MPHA.gt.60    <- which(rawData2$mass.liq.frac[ind.dmax.gt00.FZDZ.gt.00] < 0.33 & rawData2$MIXPHA.int.match.matrix.max[ind.dmax.gt00.FZDZ.gt.00] >= 0.60)
ind.mlf.lt.33.MPHA.lt.60    <- which(rawData2$mass.liq.frac[ind.dmax.gt00.FZDZ.gt.00] < 0.33 & rawData2$MIXPHA.int.match.matrix.max[ind.dmax.gt00.FZDZ.gt.00] < 0.60)
frac.mlf.gt.66.MPHA.lt.65   <- length(ind.mlf.gt.66.MPHA.lt.65) / (length(ind.mlf.gt.66.MPHA.lt.65) + length(ind.mlf.gt.66.MPHA.gt.65))
frac.mlf.lt.33.MPHA.gt.60   <- length(ind.mlf.lt.33.MPHA.gt.60) / (length(ind.mlf.lt.33.MPHA.gt.60) + length(ind.mlf.lt.33.MPHA.lt.60))
text(0.15, 1.02, paste(round(frac.mlf.lt.33.MPHA.gt.60, digits=2)*100, "% pts MLF<0.33 over 0.60 int", sep=""))
text(0.85, 0.30, paste(round(frac.mlf.gt.66.MPHA.lt.65, digits=2)*100, "% pts MLF>0.66 under 0.65 int", sep=""))

ind.mlf.gt.66               <- which(rawData2$mass.liq.frac[ind.dmax.gt00.FZDZ.gt.00] > 0.66)
ind.mlf.bt.33.66            <- which(rawData2$mass.liq.frac[ind.dmax.gt00.FZDZ.gt.00] >= 0.33 & rawData2$mass.liq.frac[ind.dmax.gt00.FZDZ.gt.00] <= 0.66)
ind.mlf.lt.33               <- which(rawData2$mass.liq.frac[ind.dmax.gt00.FZDZ.gt.00] < 0.33)
LWC.ros.median.mlf.gt.66    <- round(median(rawData2$LWC.ros[ind.dmax.gt00.FZDZ.gt.00[ind.mlf.gt.66]], na.rm=TRUE),    digits=2)
LWC.ros.median.mlf.bt.33.66 <- round(median(rawData2$LWC.ros[ind.dmax.gt00.FZDZ.gt.00[ind.mlf.bt.33.66]], na.rm=TRUE), digits=2)
LWC.ros.median.mlf.lt.33    <- round(median(rawData2$LWC.ros[ind.dmax.gt00.FZDZ.gt.00[ind.mlf.lt.33]], na.rm=TRUE),    digits=2)
LWC.ros.mean.mlf.gt.66      <- round(mean(rawData2$LWC.ros[ind.dmax.gt00.FZDZ.gt.00[ind.mlf.gt.66]], na.rm=TRUE),      digits=2)
LWC.ros.mean.mlf.bt.33.66   <- round(mean(rawData2$LWC.ros[ind.dmax.gt00.FZDZ.gt.00[ind.mlf.bt.33.66]], na.rm=TRUE),   digits=2)
LWC.ros.mean.mlf.lt.33      <- round(mean(rawData2$LWC.ros[ind.dmax.gt00.FZDZ.gt.00[ind.mlf.lt.33]], na.rm=TRUE),      digits=2)
text(0.16, 0.05, paste("LWC-median=", LWC.ros.median.mlf.lt.33,    sep=""))
text(0.50, 0.05, paste("LWC-median=", LWC.ros.median.mlf.bt.33.66, sep=""))
text(0.84, 0.05, paste("LWC-median=", LWC.ros.median.mlf.gt.66,    sep=""))

#-----------
# plot #6, column 1 : dmax versus LWC and temp versus LWC, colored by LWC-Ros
plot.new()
par(mfrow=c(1,2))
par(mar=c(11, 4, 4, 4))
plot(log(rawData2$Dmax.um[ind.cdp.gt01]), pmax(rawData2$LWC.ros[ind.cdp.gt01], rawData2$LWC.nev[ind.cdp.gt01]), pch=20, col="black", xaxt="n", xlim=c(3, 8), ylim=c(0, max(rawData2$LWC.ros, na.rm=TRUE)), main=paste(radia.date, sep=""), xlab="D-max [um]", ylab="max(LWC-Rose/LWC-CDP) [gm-3]")
axis(1, at=3:8, labels=c(20, 55, 148, 403, 1096, 2980))
text(5.7, 0.90, paste("Data represents ", as.character(num.total), " 30-s averages", sep=""))
text(4.7, 0.85, paste("over ", as.character(num.cases+2), " case dates", sep=""))
grid()

# plot #6, column 2: Temp versus LWC
plot(pmin(rawData2$T.rev.C[ind.cdp.gt01], rawData2$T.ros.C[ind.cdp.gt01]), pmax(rawData2$LWC.ros[ind.cdp.gt01], rawData2$LWC.nev[ind.cdp.gt01]), pch=20, col="black", xlim=c(-35, 2), ylim=c(0, max(rawData2$LWC.ros, na.rm=TRUE)), xlab="Temp [C]", ylab="max(LWC-Rose/LWC-Nev) [gm-3]")
grid()

# stats on fit
summary(fit)
#summary(fit2)
summary(fit3)

#----------
# plot #7
if (radia.date == radia.date) {
  ind.timeseries1                               <- which(rawData2$radia.date == radia.date & rawData2$radia.hh.num == 14 & rawData2$radia.mm.num > 3)
  ind.timeseries2                               <- which(rawData2$radia.date == radia.date & rawData2$radia.hh.num >= 15)
  #ind.timeseries3                               <- which(rawData2$radia.date == radia.date & rawData2$radia.hh.num ==22 & rawData2$radia.mm.num < 36)
  ind.timeseries                                <- as.numeric(c(ind.timeseries1, ind.timeseries2))
  #ind.timeseries                                <- as.numeric(c(ind.timeseries1, ind.timeseries2, ind.timeseries3))
  radia.hhmm                                    <- paste(as.character(rawData2$radia.hh.num[ind.timeseries]), ":", as.character(rawData2$radia.mm.num[ind.timeseries]), sep="")
  FZDZ.in.timeseries                            <- rawData2$FZDZ.int.match.matrix.max[ind.timeseries]
  #FZDZ.in.timeseries[is.na(FZDZ.in.timeseries)] <- -0.1
  MPHA.in.timeseries                            <- rawData2$MIXPHA.int.match.matrix.max[ind.timeseries]
  #MPHA.in.timeseries[MPHA.in.timeseries<=0.33]  <- NA
  
  # top row
  plot.new()
  par(mfrow=c(3,1), mar=c(5,5,2,5))
  plot(seq(1, length(rawData2$FZDZ.int.match.matrix.max[ind.timeseries])), FZDZ.in.timeseries, col="blueviolet", type="l", ylim=c(-0.1, 1.0), main=paste(radia.date, " ", radia.hhmm[1], " to ", radia.hhmm[length(ind.timeseries)], " UTC", sep=""), xlab=NA, ylab=NA)
  lines(seq(1, length(rawData2$FZDZ.int.match.matrix.max[ind.timeseries])), MPHA.in.timeseries, col="orange", type="l", xaxt="none")
  abline(h=0.6, lty=2)
  grid()
  #axis(side=1, rep(20,length(ind.timeseries)), las=2, cex.axis=0.8, font=5, col.axis="black")
  mtext(side=2, line=3, 'RadIA int')
  
  # middle row
  plot(seq(1, length(rawData2$FZDZ.int.match.matrix.max[ind.timeseries])), log(rawData2$dmax.um[ind.timeseries]), col="green", type="l", ylim=c(4, 8), axes=TRUE, xlab=NA, ylab=NA, yaxt="none")
  axis(side=2, at=c(4, 5, 6, 7, 8), labels=c(55, 148, 403, 1096, 2980), las=2, cex.axis=0.8, font=4, col.axis="green")
  mtext(side=2, line=3, 'D-99% [um]', col="green")
  grid()
  par(new=TRUE)
  plot(seq(1, length(rawData2$FZDZ.int.match.matrix.max[ind.timeseries])), rawData2$LWC.nev[ind.timeseries], axes=FALSE, col="red", type="l", lty=3, ylim=c(0.0, 0.5), xlab=NA, ylab=NA, xaxt="none")
  lines(seq(1, length(rawData2$FZDZ.int.match.matrix.max[ind.timeseries])), rawData2$LWC.ros[ind.timeseries], axes=FALSE, col="orange", type="l", lty=3)
  axis(side=4, seq(0, 0.5, 0.1),   las=2, cex.axis=0.8, font=4, col.axis="red")
  mtext(side=4, text='LWC max(NEV/CDP) [gm-3]', line=2.5, col="red")
  grid()
 
  # bottom row
  plot(seq(1, length(rawData2$FZDZ.int.match.matrix.max[ind.timeseries])), rawData2$AC.alt.m[ind.timeseries], col="black", type="l", lty=2, ylim=c(4000, 5500), xlab=NA, ylab=NA)
  mtext(side=2, line=3, 'Alt [m]')
  grid()
  par(new=TRUE)
  plot(seq(1, length(rawData2$FZDZ.int.match.matrix.max[ind.timeseries])), rawData2$T.rev.C[ind.timeseries], axes=FALSE, col="blue", type="l", lty=2, ylim=c(-26, 2), xlab=NA, ylab=NA, xaxt="none")
  axis(side=4, seq(-26, 2, 4),   las=2, cex.axis=0.8, font=4, col.axis="blue")
  mtext(side=4, text='Temp [deg C]', line=2.5, col="blue")
  grid()
  
}

# 
hist(rawData2$Dmax.um[ind.cdp.gt01])
grid()

#---------------------------------------
# individual case period analyses
#---------------------------------------

# 20170108 3:00 to 3:15 UTC
par(mfrow=c(2,1))
par(mar=c(4,4,4,4))
#ind.segment.20170108.02Z     <- which(rawData2$radia.date == "20170108" & rawData2$radia.hh.num == 2 & rawData2$radia.mm.num >= 45 & rawData2$radia.mm.num < 60)
ind.segment.20170108.03Z     <- which(rawData2$radia.date == "20170108" & rawData2$radia.hh.num == 3 & rawData2$radia.mm.num >= 0 & rawData2$radia.mm.num <= 15)
#ind.segment.20170108.02to03Z <- union(ind.segment.20170108.02Z, ind.segment.20170108.03Z)
plot(rawData2$LWC.ros[ind.segment.20170108.03Z],                      col="red",    pch=20, type="b", ylim=c(0,1), main="20170108 03:00 to 03:15 UTC", xlab="Time [UTC]", ylab="LWC/INT/FRAC")
lines(rawData2$LWC.nev[ind.segment.20170108.02to03Z],                 col="red",    pch=21, type="b")
lines(rawData2$SLW.int.match.matrix.max[ind.segment.20170108.03Z],    col="green",  pch=21, type="l")
lines(rawData2$FZDZ.int.match.matrix.max[ind.segment.20170108.03Z],   col="blue",   pch=21, type="l")
lines(rawData2$MIXPHA.int.match.matrix.max[ind.segment.20170108.03Z], col="orange", pch=21, type="l")
lines(1-rawData2$mass.liq.frac[ind.segment.20170108.03Z],             col="grey",   pch=20, type="l", lwd=3)
abline(h=0.59, lty=2)
grid()
legend(24, 0.8, legend=c("LWC-Ros", "LWC-Nev", "MIF", "small-drop", "large-drop", "mixed-phase"), col=c("red", "red", "grey", "green", "blue", "orange"), pch=c(20, 21, 20, 21, 21, 21))
plot(rawData2$T.ros.C[ind.segment.20170108.03Z], col="red", pch=20, type="b", ylim=c(-20,0), xlab="Time [UTC]", ylab="Temp [deg C]")
grid()

# 20170109 04:36 to 04:46 UTC
par(mfrow=c(2,1))
par(mar=c(4,4,4,4))
ind.segment.20170109.04Z   <- which(rawData2$radia.date == "20170109" & rawData2$radia.hh.num == 4 & rawData2$radia.mm.num >= 32 & rawData2$radia.mm.num <= 44)
plot(rawData2$LWC.ros[ind.segment.20170109.04Z[1:25]],                      col="red",    pch=20, type="b", ylim=c(0,1), main="20170109 04:32 to 04:44 UTC", xlab="Time [UTC]", ylab="LWC/INT/FRAC")
lines(rawData2$LWC.nev[ind.segment.20170109.04Z[1:25]]-0.02,                col="red",    pch=21, type="b")
lines(rawData2$SLW.int.match.column.max[ind.segment.20170109.04Z[1:25]],    col="green",  pch=21, type="l")
lines(rawData2$FZDZ.int.match.column.max[ind.segment.20170109.04Z[1:25]],   col="blue",   pch=21, type="l")
lines(rawData2$MIXPHA.int.match.column.max[ind.segment.20170109.04Z[1:25]], col="orange", pch=21, type="l")
lines(1-rawData2$mass.liq.frac[ind.segment.20170109.04Z[1:25]],             col="grey",   pch=21, type="l", lwd=3)
abline(h=0.60, lty=2)
grid()
legend(20, 0.55, legend=c("LWC-Ros", "LWC-Nev", "MIF", "small-drop", "large-drop", "mixed-phase"), col=c("red", "red", "grey", "green", "blue", "orange"), pch=c(20,21,20,21,21,21))
plot(rawData2$T.ros.C[ind.segment.20170109.04Z[1:25]], col="red", pch=20, type="b", ylim=c(-20,0), xlab="Time [UTC]", ylab="Temp [deg C]")
grid()

# 20170109 6:31 to 6:50 UTC very low LWCs
par(mfrow=c(2,1))
par(mar=c(4,4,4,4))
ind.segment.20170109.06Z   <- which(rawData2$radia.date == "20170109" & rawData2$radia.hh.num == 6 & rawData2$radia.mm.num >= 31 & rawData2$radia.mm.num <= 50)
plot(rawData2$LWC.ros[ind.segment.20170109.06Z], col="red", pch=20, ylim=c(0,1), main="20170109 06:31 to 06:50 UTC", xlab="Time [UTC]", ylab="LWC/INT/FRAC")
lines(rawData2$LWC.nev[ind.segment.20170109.06Z, col="red", pch=21])
lines(rawData2$SLW.int.match.column.max[ind.segment.20170109.06Z],    col="green",  pch=21, type="l")
lines(rawData2$FZDZ.int.match.column.max[ind.segment.20170109.06Z],   col="blue",   pch=21, type="l")
lines(rawData2$MIXPHA.int.match.column.max[ind.segment.20170109.06Z], col="orange", pch=21, type="l")
lines(1-rawData2$mass.liq.frac[ind.segment.20170109.06Z],             col="grey",   pch=21, type="l", lwd=3)
abline(h=0.60, lty=2)
grid()
legend()
plot(rawData2$T.ros.C[ind.segment.20170109.06Z], col="red", pch=20, ylim=c(-20,0), xlab="Time [UTC]", ylab="Temp [deg C]")
grid()

# 20170122 21:50 to 22:24 UTC
par(mfrow=c(2,1))
par(mar=c(4,4,4,4))
ind.segment.20170122.21Z     <- which(rawData2$radia.date == "20170122" & rawData2$radia.hh.num == 21 & rawData2$radia.mm.num > 50)
ind.segment.20170122.22Z     <- which(rawData2$radia.date == "20170122" & rawData2$radia.hh.num == 22 & rawData2$radia.mm.num <= 24)
ind.segment.20170122.21to22Z <- union(ind.segment.20170122.21Z, ind.segment.20170122.22Z)
plot(rawData2$LWC.ros[ind.segment.20170122.21to22Z],                      col="red",    pch=20, type="b", ylim=c(0,1), main="20170122 21:50 to 22:24 UTC", xlab="Time [UTC]", ylab="LWC/INT/FRAC")
lines(rawData2$LWC.nev[ind.segment.20170122.21to22Z],                     col="red",    pch=21, type="b")
lines(rawData2$SLW.int.match.matrix.max[ind.segment.20170122.21to22Z],    col="green",  pch=21, type="l")
lines(rawData2$FZDZ.int.match.matrix.max[ind.segment.20170122.21to22Z],   col="blue",   pch=21, type="l")
lines(rawData2$MIXPHA.int.match.matrix.max[ind.segment.20170122.21to22Z], col="orange", pch=21, type="l")
lines(1-rawData2$mass.liq.frac[ind.segment.20170122.21to22Z],             col="grey",   pch=20, type="l", lwd=3)
grid()
abline(h=0.6, lty=2)
legend(55, 0.7, legend=c("LWC-Ros", "LWC-Nev", "MIF", "small-drop", "large drop", "mixed-phase"), col=c("red", "red", "grey", "green", "blue", "orange"), pch=c(20,21,20,21,21,21))
plot(rawData2$T.ros.C[ind.segment.20170122.21to22Z], col="red", pch=20, type="b", ylim=c(-15, 0), xlab="Time [UTC]", ylab="Temp [deg C]")
grid()

# 20170131 20:21 to 20:25 UTC
par(mfrow=c(2,1))
par(mar=c(4,4,4,4))
ind.segment.20170131.20Z   <- which(rawData2$radia.date == "20170131" & rawData2$radia.hh.num == 20 & rawData2$radia.mm.num >= 19 & rawData2$radia.mm.num <= 24)
plot(rawData2$LWC.ros[ind.segment.20170131.20Z],                      col="red",    pch=20, type="b", ylim=c(0,1), main="20170131 20:19 to 20:24 UTC", xlab="time [UTC]", ylab="LWC/INT/FRAC")
lines(rawData2$LWC.nev[ind.segment.20170131.20Z],                     col="red",    pch=21, type="b")
lines(rawData2$SLW.int.match.column.max[ind.segment.20170131.20Z],    col="green",  pch=21, type="l")
lines(rawData2$FZDZ.int.match.column.max[ind.segment.20170131.20Z],   col="blue",   pch=21, type="l")
lines(rawData2$MIXPHA.int.match.column.max[ind.segment.20170131.20Z], col="orange", pch=21, type="l")
lines(1-rawData2$mass.liq.frac[ind.segment.20170131.20Z],             col="grey",   pch=20, type="l", lwd=3)
grid()
abline(h=0.59, lty=2)
legend(10.5, 1.0, legend=c("LWC-Ros", "LWC-Nev", "MIF", "small-drop", "large-drop", "mixed-phase"), col=c("red", "red", "grey", "green", "blue", "orange"), pch=c(20,21,20,21,21,21))
plot(rawData2$T.ros.C[ind.segment.20170131.20Z], col="red", pch=20, type="b", ylim=c(-15,0), xlab="Time [UTC]", ylab="Temp [deg C]")
grid()

# 20170131 21:15 to 21:23 UTC
par(mfrow=c(2,1))
par(mar=c(4,4,4,4))
ind.segment.20170131.21Z   <- which(rawData2$radia.date == "20170131" & rawData2$radia.hh.num == 21 & rawData2$radia.mm.num >= 15 & rawData2$radia.mm.num <= 23)
plot(rawData2$LWC.ros[ind.segment.20170131.21Z],                      col="red",    pch=20, type="b", ylim=c(0,1), main="20170131 21:15 to 21:23 UTC", xlab="Time [UTC]", ylab="LWC/INT/FRAC")
lines(rawData2$LWC.nev[ind.segment.20170131.21Z],                     col="red",    pch=21, type="b")
lines(rawData2$SLW.int.match.column.max[ind.segment.20170131.21Z],    col="green",  pch=21, type="l")
lines(rawData2$FZDZ.int.match.column.max[ind.segment.20170131.21Z],   col="blue",   pch=21, type="l")
lines(rawData2$MIXPHA.int.match.column.max[ind.segment.20170131.21Z], col="orange", pch=21, type="l")
lines(1-rawData2$mass.liq.frac[ind.segment.20170131.21Z], col="grey", pch=20, type="l", lwd=3)
grid()
abline(h=0.6, lty=2)
legend(16.5, 1.0, legend=c("LWC-Ros", "LWC-Nev", "MIF", "small-drop", "large-drop", "mixed-phase"), col=c("red", "red", "grey", "green", "blue", "orange"), pch=c(20,21,20,21,21,21))
plot(rawData2$T.ros.C[ind.segment.20170131.21Z], col="red", pch=20, type="b", ylim=c(-15,0), xlab="Time [UTC]", ylab="Temp [deg C]")
grid()

# plot massliqfrac and LWCROS versus FZDZ and MIXPHA for 1/22/2017
par(mfrow=c(2, 1))
par(mar=c(4, 4, 4, 4))
plot(rawData2$mass.liq.frac[ind.dmax.gt00.FZDZ.gt.00[219:338]], type="o", lty=1, pch=20, main="20170122")
lines(rawData2$LWC.ros[ind.dmax.gt00.FZDZ.gt.00[219:338]], type="o", col="red")
lines(rawData2$MIXPHA.int.match.matrix.max[ind.dmax.gt00.FZDZ.gt.00[219:338]], type="o", col="blue")
lines(rawData2$FZDZ.int.match.matrix.max[ind.dmax.gt00.FZDZ.gt.00[219:338]], type="o", col="green")
grid()
plot(rawData2$dmax.um[ind.dmax.gt00.FZDZ.gt.00[219:338]], type="o", lty=1, pch=20, ylim=c(0,500))
grid()

# plot massliqfrac and LWCROS versus FZDZ and MIXPHA for 1/31/2017
par(mfrow=c(2, 1))
par(mar=c(4, 4, 4, 4))
plot(rawData2$mass.liq.frac[ind.dmax.gt00.FZDZ.gt.00[339:474]], type="o", lty=1, pch=20, main="20170131")
lines(rawData2$LWC.ros[ind.dmax.gt00.FZDZ.gt.00[339:474]], type="l", col="red")
lines(rawData2$LWC.nev[ind.dmax.gt00.FZDZ.gt.00[339:474]], type="l", col="red", lty=2)
lines(rawData2$MIXPHA.int.match.matrix.max[ind.dmax.gt00.FZDZ.gt.00[339:474]], type="o", col="blue")
lines(rawData2$FZDZ.int.match.matrix.max[ind.dmax.gt00.FZDZ.gt.00[339:474]], type="o", col="green")
grid()
plot(rawData2$dmax.um[ind.dmax.gt00.FZDZ.gt.00[339:474]], type="o", lty=1, pch=20, ylim=c(0,100))
grid()

#############################
# PERIOD FOR RADAR CONF PAPER
#############################
# plot massliqfrac and LWCROS versus FZDZ and MIXPHA for 1/31/2017 inds #1-12 (20:19-20:25Z)
par(mfrow=c(2, 1))
par(mar=c(4, 4, 4, 4))
i.srt <- 148
i.stp <- 255
#i.stp <- length(ind.dmax.gt00.FZDZ.gt.00)
plot(1-rawData2$mass.liq.frac[i.srt:i.stp],                type="o", col="grey", lty=1, pch=20, main=paste("20170131 ", rawData2$radia.hh.num[i.srt], ":", rawData2$radia.mm.num[i.srt], "-",  rawData2$radia.hh.num[i.stp], ":", rawData2$radia.mm.num[i.stp], " UTC", sep=""), ylim=c(0,1), xlab="", ylab="")
lines(rawData2$LWC.ros[i.srt:i.stp],                     type="l", col="red", lty=1)
lines(rawData2$LWC.nev[i.srt:i.stp],                     type="l", col="red", lty=2)
#lines(rawData2$LWC.cdp[ind.dmax.gt00.FZDZ.gt.00[i.srt:i.stp]]-0.07,                type="l", col="red", lty=3)
lines(rawData2$MIXPHA.int.match.matrix.max[i.srt:i.stp], type="l", col="orange")
lines(rawData2$SLW.int.match.matrix.mean[i.srt:i.stp],   type="l", col="blue")
lines(rawData2$FZDZ.int.match.matrix.max[i.srt:i.stp],   type="l", col="green")
par(xpd=TRUE)
legend(5.0,  -0.06, legend=c("MASS-LQ-FR", "FZDZ int", "SLW-int", "MPHA-int"), col=c("black", "green", "blue", "orange"), lty=c(1, 1, 1, 1), cex=0.8)
legend(35.0, -0.06, legend=c("LWC-ROSE", "LWC-NEV"), col=c("red", "red"), lty=c(1, 2), cex=0.8)
x.ticks <- rawData2$radia.mm.num[1:12]
x.ticks <- c(19, 19, 20, 20, 21, 21, 22, 22, 23, 23, 24, 24)
axis(side=1, at=x.ticks, labels=FALSE)
text(x=x.ticks,  par("usr")[3], labels = x.ticks, srt = 45, pos = 1, xpd = TRUE)
grid()
plot(rawData2$dmax.um[i.srt:i.stp], type="o", lty=1, pch=20, ylim=c(0,100))
grid()

# plot massliqfrac and LWCROS versus FZDZ and MIXPHA for 2/7/2017
par(mfrow=c(2, 1))
par(mar=c(4, 4, 4, 4))
plot(rawData2$mass.liq.frac[ind.dmax.gt00.FZDZ.gt.00[494:540]], type="o", lty=1, pch=20, main="20170207")
lines(rawData2$LWC.ros[ind.dmax.gt00.FZDZ.gt.00[494:540]], type="o", col="red")
lines(rawData2$MIXPHA.int.match.matrix.max[ind.dmax.gt00.FZDZ.gt.00[494:540]], type="o", col="blue")
lines(rawData2$FZDZ.int.match.matrix.max[ind.dmax.gt00.FZDZ.gt.00[494:540]], type="o", col="green")
grid()
plot(rawData2$dmax.um[ind.dmax.gt00.FZDZ.gt.00[494:540]], type="o", lty=1, pch=20, ylim=c(0,100))
grid()

# plot massliqfrac and LWCROS versus FZDZ and MIXPHA for 2/21/2017
par(mfrow=c(2, 1))
par(mar=c(4, 4, 4, 4))
plot(rawData2$mass.liq.frac[ind.dmax.gt00.FZDZ.gt.00[546:641]], type="o", lty=1, pch=20, main="20170221")
lines(rawData2$LWC.ros[ind.dmax.gt00.FZDZ.gt.00[546:641]], type="o", col="red")
lines(rawData2$MIXPHA.int.match.matrix.max[ind.dmax.gt00.FZDZ.gt.00[546:641]], type="o", col="blue")
lines(rawData2$FZDZ.int.match.matrix.max[ind.dmax.gt00.FZDZ.gt.00[546:641]], type="o", col="green")
grid()
plot(rawData2$dmax.um[ind.dmax.gt00.FZDZ.gt.00[546:641]], type="o", lty=1, pch=20, ylim=c(0,100))
grid()

# plot massliqfrac and LWCROS versus FZDZ and MIXPHA for 3/9/2017
par(mfrow=c(2, 1))
par(mar=c(4, 4, 4, 4))
plot(rawData2$mass.liq.frac[ind.dmax.gt00.FZDZ.gt.00[642:841]], type="o", lty=1, pch=20, main="20170309")
lines(rawData2$LWC.ros[ind.dmax.gt00.FZDZ.gt.00[642:841]], type="o", col="red")
lines(rawData2$MIXPHA.int.match.matrix.max[ind.dmax.gt00.FZDZ.gt.00[642:841]], type="o", col="blue")
lines(rawData2$FZDZ.int.match.matrix.max[ind.dmax.gt00.FZDZ.gt.00[642:841]], type="o", col="green")
grid()
plot(rawData2$dmax.um[ind.dmax.gt00.FZDZ.gt.00[642:841]], type="o", lty=1, pch=20, ylim=c(0,500))
grid()

# plot massliqfrac and LWCROS versus FZDZ and MIXPHA for 3/9/2017 first 13 inds
par(mfrow=c(2, 1))
par(mar=c(4, 4, 4, 4))
plot(rawData2$mass.liq.frac[ind.dmax.gt00.FZDZ.gt.00[642:654]], type="o", lty=1, pch=20, main="20170309", ylim=c(0, 1))
lines(rawData2$LWC.ros[ind.dmax.gt00.FZDZ.gt.00[642:654]], type="o", col="red")
lines(rawData2$LWC.nev[ind.dmax.gt00.FZDZ.gt.00[642:654]], type="o", col="red", lty=2)
lines(rawData2$LWC.cdp[ind.dmax.gt00.FZDZ.gt.00[642:654]], type="o", col="red", lty=3)
lines(rawData2$MIXPHA.int.match.matrix.max[ind.dmax.gt00.FZDZ.gt.00[642:654]], type="o", col="blue")
lines(rawData2$FZDZ.int.match.matrix.max[ind.dmax.gt00.FZDZ.gt.00[642:654]], type="o", col="green")
grid()
plot(rawData2$dmax.um[ind.dmax.gt00.FZDZ.gt.00[642:654]], type="o", lty=1, pch=20, ylim=c(0,110))
grid()

# terrain blockage plots
# terrain height above radar within 25 km vs PODY
bin.height.terrain.above.radar <- c(62.50, 182.50, 307.50, 432.50, 557.50)
bin.PODY                       <- c( 0.72,   0.83,   0.74,   0.45,   0.15  )
plot(bin.height.terrain.above.radar, bin.PODY, type="o", lty=1, pch=20, ylim=c(0,1), xlab="terrain height above radar within 25km [m]", ylab="RadIA POD-Y")
grid()

#
par(mfrow=c(1, 1))
par(mar=c(4, 4, 4, 4))
#terr.alt.NE     <- c(1050,  1150,    1200,   1250,  1350,  1280,  1300,  1400,  1380, 1370)
#azim.NE         <- c(  0,      5,      10,     15,    20,    25,    30,    35,    40,   45)
terr.alt.NW     <- c(940,    960,     940,    960,   980,   990,  1100,  1100,  1050)
azim.NW         <- c(320,    325,     330,    335,   340,   345,   350,   355,   360)
#ind.NE.FZDZge00 <- which(rawData2$round.azim.radar.to.ac.deg..digits...2.[ind.dmax.gt100.cdp.gt01.mlf.gt20] <= 50 & rawData2$FZDZ.int.match.matrix.max[ind.dmax.gt100.cdp.gt01.mlf.gt20] >= 0)
ind.NW.FZDZge00 <- which(rawData2$round.azim.radar.to.ac.deg..digits...2.[ind.dmax.gt100.cdp.gt01.mlf.gt20] > 320 & rawData2$FZDZ.int.match.matrix.max[ind.dmax.gt100.cdp.gt01.mlf.gt20] >= 0)
#ind.NE.FZDZna   <- which(rawData2$round.azim.radar.to.ac.deg..digits...2.[ind.dmax.gt100.cdp.gt01.mlf.gt20] <= 50 & is.na(rawData2$FZDZ.int.match.matrix.max[ind.dmax.gt100.cdp.gt01.mlf.gt20]))
ind.NW.FZDZna   <- which(rawData2$round.azim.radar.to.ac.deg..digits...2.[ind.dmax.gt100.cdp.gt01.mlf.gt20] > 320 & is.na(rawData2$FZDZ.int.match.matrix.max[ind.dmax.gt100.cdp.gt01.mlf.gt20]))
#plot(rawData2$round.azim.radar.to.ac.deg..digits...2.[ind.dmax.gt100.cdp.gt01.mlf.gt20[ind.NE.FZDZge00]], rawData2$AC.alt.m[ind.dmax.gt100.cdp.gt01.mlf.gt20[ind.NE.FZDZge00]] , col="black", type="p", lty=1, pch=21, ylim=c(0, max(rawData2$AC.alt.m[ind.ac.appO])), xlab="AZIM [deg] from KCBX", ylab="Altitude [m, MSL]")
#lines(rawData2$round.azim.radar.to.ac.deg..digits...2.[ind.dmax.gt100.cdp.gt01.mlf.gt20[ind.NE.FZDZna]], rawData2$AC.alt.m[ind.dmax.gt100.cdp.gt01.mlf.gt20[ind.NE.FZDZna]], col="red", type="p", lty=1, pch=20)
#lines(azim.NE, terr.alt.NE, type="l", col="chocolate")
plot(rawData2$round.azim.radar.to.ac.deg..digits...2.[ind.dmax.gt100.cdp.gt01.mlf.gt20[ind.NW.FZDZge00]], rawData2$AC.alt.m[ind.dmax.gt100.cdp.gt01.mlf.gt20[ind.NW.FZDZge00]] , col="black", type="p", lty=1, pch=21, ylim=c(0, max(rawData2$AC.alt.m[ind.ac.appO])), xlab="AZIM [deg] from KCBX", ylab="Altitude [m, MSL]", xlim=c(320, 360))
lines(rawData2$round.azim.radar.to.ac.deg..digits...2.[ind.dmax.gt100.cdp.gt01.mlf.gt20[ind.NW.FZDZna]], rawData2$AC.alt.m[ind.dmax.gt100.cdp.gt01.mlf.gt20[ind.NW.FZDZna]], col="red", type="p", lty=1, pch=20)
lines(azim.NW, terr.alt.NW, type="l", col="chocolate")
abline(h=968, lty=2)
text(10, 1050, "KCBX elevation")
legend(30, 3000, legend=c("FZDZ >= 0", "FZDZ = NA "), col=c("black", "red"), pch=c(21, 20))
grid()
