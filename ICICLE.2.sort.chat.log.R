#---------------------------------
# Name:     ICICLE.2.sort.chat.log.R
#
# Purpose:  NCAR scientists won't have even quicklook data until at least september 2019.
#           In order to make progress with algorithm skill analysis, the following process is conducted:
#
#           Execute the following scripts in the following order:
#
#           1. ICICLE.1.flhours.vs.days.R
#               a. builds 'icicle.df' and 'icicle.chat.ificat.df' data frames
#           2. ICICLE.2.sort.chat.log.R
#           	a. Search through CNV flight chat log to find specified keyword occurrences, then find date/time of occurrences
#           	b. If keyword occurrence is within primary radar domain (within ~100km), search CNV flight track log
#                for the lat/lon/alt/tempa/tempd values.  
#               c. Plot these metadata for individual flights and whole campaign
#               d. Save off csv files of these data from every flight
#           3. ICICLE.3.load.manip.radia.int.R
#               a. Find closest radar volume in time and find the closest RadIA INTs to the keyword occurrences.
#               b. Plot statistics on num occur, alts, mean/median/boxwhisker plots of INT values, etc
#
#
#           Approx number of hours of data to assess: 27hr and 43min plus primary time from F7,17,18,23, & 26
#           Number of flights to assess:              26 research flights, 19 with processed RadIA
#
# Created:  3.7.2019 dserke
#---------------------------------
  
#---------------------------------
# define flag variables and values
output.plot.flag <- 1
output.data.flag <- 1
#---------------------------------

#---------------------------------
# define directories 
#---------------------------------
#  from where the flight csv files live in
icicle.chat.dataframecsv.dir      <- file.path("/d1/serke/projects/case_studies/ICICLE_2018/data/chat_planet/")
icicle.flightleg.dataframecsv.dir <- file.path("/d1/serke/projects/case_studies/ICICLE_2018/data/flight_track_planet/")
#  for output plots and data
output.plot.dir                   <- file.path("/d1/serke/projects/case_studies/ICICLE_2018/images/auto_sort_chat_keyword/")
output.data.dir                   <- file.path("/d1/serke/projects/case_studies/ICICLE_2018/data/auto_sort_chat_keyword/")
#---------------------------------

#---------------------------------
# build chat and flight track file names from only ICICLE flight number
#---------------------------------
ICICLE.fl.num          <- 6
ind.fl.num             <- which(icicle.table.df$fl.num == ICICLE.fl.num)
ICICLE.track.file.name <- paste("F", as.character(ICICLE.fl.num), "_NRC_track_", as.character(icicle.2.df$date.posix[ind.fl.num]), ".csv", sep="")
ICICLE.chat.file.name  <- paste("F", as.character(ICICLE.fl.num), "_NRC_chat_",  as.character(icicle.2.df$date.posix[ind.fl.num]), ".csv", sep="")
ICICLE.track.file.name
ICICLE.chat.file.name
#---------------------------------

#---------------------------------
# load chat & flight track csv format data files
#---------------------------------
#   load chat csv file
chat.data           <- read.csv(paste(icicle.chat.dataframecsv.dir, ICICLE.chat.file.name,    sep=""))
colnames(chat.data) <- c("CHAT", "POSIX.DATE.1", "POSIX.DATE.2", "PLANET", "CAMPAIGN", "CHANNEL", "USER", "COMMENT")
# NOTE: KEY TO CHANNELS- 1=FLIGHTSCIENTIST, 2=OPSDIRECTION, 3=G2G, 4=FAA, 5=NCAR, 6=RADAR

#   load flight track csv file
flighttrack         <- read.csv(paste(icicle.flightleg.dataframecsv.dir, ICICLE.track.file.name, sep=""))
colnames(flighttrack)<- c("Type", "ac.date.posix", "lat", "lon", "alt.msl.m", "alt.wgs.m", "pres.alt.ft", "radar.alt.ft", "grnd.spd.mps", "tru.spd.mps", "ind.spd.kts", "mach.num", "vert.vel.mps", "tru.hdg.deg", "track.deg", "drift.deg", "pitch.deg", "roll.deg", "side.slip.deg", "ang.atk.deg", "amb.temp.degC", "dew.temp.degC", "tot.temp.deg", "stat.pres.mb", "dyn.pres.mb", "cab.pres.mb", "wind.spd.mps", "wind.dir.mps") # from email from M. Bastian NRC, 2/4/2019
#---------------------------------

#---------------------------------
# textual keywords for later indexing
#---------------------------------
#   define keywords of interest for each type of IFI
keywords.MIXPHA     <- c("MIXPHA",    "MP",               "mixpha",           "mixed phase", "ice liquid", "mix of ice")
keywords.AppC       <- c("smalldrop", "app c",            "AppC",             "App C",       "Appx C",     "all small drop")
keywords.AppOFZRA   <- c("FZRA",      "freezing rain")
keywords.AppOFZDZ   <- c("FZDZ",      "DZ",               "freezing drizzle", "drizzle",     "drizzling",   "SLD",        "fr dz",      "fr driz")
keywords.ALLICE     <- c("glaciated", "all ice")
keywords.DENDS      <- c("DENDRITE",  "dendrite",         "dend")
keywords.PLATES     <- c("PLATE",     "plate")

#   define scientists of interest, only use comments from them
keywords.FLSCIENTIST<- c("koroleve",  "nichmanl",         "woldem",            "cuongn",      "thompsong")
#---------------------------------

#---------------------------------
# build indices
#---------------------------------
#   index for keywords of interest
#   add one to index, seems to be C based structure where first value is zeroth index
index.MIXPHA        <- grep(paste(keywords.MIXPHA,      collapse = "|"), chat.data$COMMENT) + 1 
index.AppC          <- grep(paste(keywords.AppC,        collapse = "|"), chat.data$COMMENT) + 1                        
index.AppOFZRA      <- grep(paste(keywords.AppOFZRA,    collapse = "|"), chat.data$COMMENT) + 1                         
index.AppOFZDZ      <- grep(paste(keywords.AppOFZDZ,    collapse = "|"), chat.data$COMMENT) + 1 
index.ALLICE        <- grep(paste(keywords.ALLICE,      collapse = "|"), chat.data$COMMENT) + 1 
index.DENDS         <- grep(paste(keywords.DENDS,       collapse = "|"), chat.data$COMMENT) + 1 
index.PLATES        <- grep(paste(keywords.PLATES,      collapse = "|"), chat.data$COMMENT) + 1
index.FLSCIENTIST   <- grep(paste(keywords.FLSCIENTIST, collapse = "|"), chat.data$USER)    + 1

#   intersected indices, FLSCIENTIST plus KEYWORD lists
index.MIXPHA.I1     <- intersect(index.MIXPHA,   index.FLSCIENTIST)
index.AppC.I1       <- intersect(index.AppC,     index.FLSCIENTIST)
index.AppOFZRA.I1   <- intersect(index.AppOFZRA, index.FLSCIENTIST)
index.AppOFZDZ.I1   <- intersect(index.AppOFZDZ, index.FLSCIENTIST)
index.ALLICE.I1     <- intersect(index.ALLICE,   index.FLSCIENTIST)
index.DENDS.I1      <- intersect(index.DENDS,    index.FLSCIENTIST)
index.PLATES.I1     <- intersect(index.PLATES,   index.FLSCIENTIST)
length(index.MIXPHA.I1)
length(index.AppC.I1)
length(index.AppOFZRA.I1)
length(index.AppOFZDZ.I1)
length(index.ALLICE.I1)
length(index.DENDS.I1)
length(index.PLATES.I1)

#   index of times within primary radar domain, from icicle.2.df data frame, which originates in script named 'ICICLE.flhours.vs.days.R'
index.PRIRADARDOM   <- which(as.POSIXct(chat.data$POSIX.DATE.1) >= as.POSIXct(paste(icicle.2.df$date.posix[ind.fl.num], " ", icicle.2.df$pri.radar.srt[ind.fl.num], ":00", sep='')) & as.POSIXct(chat.data$POSIX.DATE.1) <= as.POSIXct(paste(icicle.2.df$date.posix[ind.fl.num], " ", icicle.2.df$pri.radar.end[ind.fl.num], ":00", sep=""))) + 1

#   intersected indices, FLSCIENTIST/KEYWORD (equivalent to *.I1 designation) plus time in primary radar domain (index.PRIRADARDOM)
index.MIXPHA.I2     <- intersect(index.MIXPHA.I1,   index.PRIRADARDOM)
index.AppC.I2       <- intersect(index.AppC.I1,     index.PRIRADARDOM)
index.AppOFZRA.I2   <- intersect(index.AppOFZRA.I1, index.PRIRADARDOM)
index.AppOFZDZ.I2   <- intersect(index.AppOFZDZ.I1, index.PRIRADARDOM)
index.ALLICE.I2     <- intersect(index.ALLICE.I1,   index.PRIRADARDOM)
index.DENDS.I2      <- intersect(index.DENDS.I1,    index.PRIRADARDOM)
index.PLATES.I2     <- intersect(index.PLATES.I1,   index.PRIRADARDOM)

#   index the keyword occurrence back to time
time.MIXPHA         <- chat.data$POSIX.DATE.1[index.MIXPHA.I2]
time.AppC           <- chat.data$POSIX.DATE.1[index.AppC.I2]
time.AppOFZRA       <- chat.data$POSIX.DATE.1[index.AppOFZRA.I2]
time.AppOFZDZ       <- chat.data$POSIX.DATE.1[index.AppOFZDZ.I2]
time.ALLICE         <- chat.data$POSIX.DATE.1[index.ALLICE.I2]
time.DENDS          <- chat.data$POSIX.DATE.1[index.DENDS.I2]
time.PLATES         <- chat.data$POSIX.DATE.1[index.PLATES.I2]
length(time.MIXPHA)
length(time.AppC)
length(time.AppOFZRA)
length(time.AppOFZDZ)
length(time.ALLICE)
length(time.DENDS)
length(time.PLATES)

#-------------------------------------------------

#-------------------------------------------------
# loop through all times when MIXPHA keyword was called out by flight scientist
if (length(time.MIXPHA) > 0) {
  ind.FLTRACKTIME.4MIXPHA <- rep(0, length(time.MIXPHA))
  time.FLTRACK.4MIXPHA    <- rep(0, length(time.MIXPHA))
  lat.FLTRACK.4MIXPHA     <- rep(0, length(time.MIXPHA))
  lon.FLTRACK.4MIXPHA     <- rep(0, length(time.MIXPHA))
  alt.FLTRACK.4MIXPHA     <- rep(0, length(time.MIXPHA))
  tempa.FLTRACK.4MIXPHA   <- rep(0, length(time.MIXPHA))
  tempd.FLTRACK.4MIXPHA   <- rep(0, length(time.MIXPHA))
  comnt.FLTRACK.4MIXPHA   <- rep(0, length(time.MIXPHA))
  for (i in 1:length(time.MIXPHA)) {
    # find index of when chat time matches flight track time
    ind.FLTRACKTIME.4MIXPHA[i] <- which(abs(difftime(flighttrack$ac.date.posix, time.MIXPHA[i])) == min(abs(difftime(flighttrack$ac.date.posix, time.MIXPHA[i])))) + 1
    # record time/lat/lon/alt/temp at AC location for time when IFI keyword was observed by flight scientist
    time.FLTRACK.4MIXPHA[i]    <- flighttrack$ac.date.posix[ind.FLTRACKTIME.4MIXPHA[i]]
    lat.FLTRACK.4MIXPHA[i]     <- flighttrack$lat[ind.FLTRACKTIME.4MIXPHA[i]]
    lon.FLTRACK.4MIXPHA[i]     <- flighttrack$lon[ind.FLTRACKTIME.4MIXPHA[i]]
    alt.FLTRACK.4MIXPHA[i]     <- flighttrack$pres.alt.ft[ind.FLTRACKTIME.4MIXPHA[i]]
    tempa.FLTRACK.4MIXPHA[i]   <- flighttrack$amb.temp.degC[ind.FLTRACKTIME.4MIXPHA[i]]
    tempd.FLTRACK.4MIXPHA[i]   <- flighttrack$dew.temp.degC[ind.FLTRACKTIME.4MIXPHA[i]]
    comnt.FLTRACK.4MIXPHA[i]   <- as.character(chat.data$COMMENT[index.MIXPHA.I2[i] - 1])
  }
} else {
  time.FLTRACK.4MIXPHA    <- -99
  lat.FLTRACK.4MIXPHA     <- -99
  lon.FLTRACK.4MIXPHA     <- -99
  alt.FLTRACK.4MIXPHA     <- -99
  tempa.FLTRACK.4MIXPHA   <- -99
  tempd.FLTRACK.4MIXPHA   <- -99
  comnt.FLTRACK.4MIXPHA   <- -99
}

# loop through all times when AppC keyword was called out by flight scientist
if (length(time.AppC) > 0) {
  ind.FLTRACKTIME.4AppC <- rep(0, length(time.AppC))
  time.FLTRACK.4AppC    <- rep(0, length(time.AppC))
  lat.FLTRACK.4AppC     <- rep(0, length(time.AppC))
  lon.FLTRACK.4AppC     <- rep(0, length(time.AppC))
  alt.FLTRACK.4AppC     <- rep(0, length(time.AppC))
  tempa.FLTRACK.4AppC   <- rep(0, length(time.AppC))
  tempd.FLTRACK.4AppC   <- rep(0, length(time.AppC))
  comnt.FLTRACK.4AppC   <- rep(0, length(time.AppC))
  for (i in 1:length(time.AppC)) {
    # find index of when chat time matches flight track time
    ind.FLTRACKTIME.4AppC[i] <- which(abs(difftime(flighttrack$ac.date.posix, time.AppC[i])) == min(abs(difftime(flighttrack$ac.date.posix, time.AppC[i])))) + 1
    # record time/lat/lon/alt/temp at AC location for time when IFI keyword was observed by flight scientist
    time.FLTRACK.4AppC[i]    <- flighttrack$ac.date.posix[ind.FLTRACKTIME.4AppC[i]]
    lat.FLTRACK.4AppC[i]     <- flighttrack$lat[ind.FLTRACKTIME.4AppC[i]]
    lon.FLTRACK.4AppC[i]     <- flighttrack$lon[ind.FLTRACKTIME.4AppC[i]]
    alt.FLTRACK.4AppC[i]     <- flighttrack$pres.alt.ft[ind.FLTRACKTIME.4AppC[i]]
    tempa.FLTRACK.4AppC[i]   <- flighttrack$amb.temp.degC[ind.FLTRACKTIME.4AppC[i]]
    tempd.FLTRACK.4AppC[i]   <- flighttrack$dew.temp.degC[ind.FLTRACKTIME.4AppC[i]]
    comnt.FLTRACK.4AppC[i]   <- as.character(chat.data$COMMENT[index.AppC.I2[i] - 1])
  }
}  else {
  time.FLTRACK.4AppC    <- -99
  lat.FLTRACK.4AppC     <- -99
  lon.FLTRACK.4AppC     <- -99
  alt.FLTRACK.4AppC     <- -99
  tempa.FLTRACK.4AppC   <- -99
  tempd.FLTRACK.4AppC   <- -99
  comnt.FLTRACK.4AppC   <- -99
}

# loop through all times when AppOFZRA keyword was called out by flight scientist
if (length(time.AppOFZRA) > 0) {
  ind.FLTRACKTIME.4AppOFZRA <- rep(0, length(time.AppOFZRA))
  time.FLTRACK.4AppOFZRA    <- rep(0, length(time.AppOFZRA))
  lat.FLTRACK.4AppOFZRA     <- rep(0, length(time.AppOFZRA))
  lon.FLTRACK.4AppOFZRA     <- rep(0, length(time.AppOFZRA))
  alt.FLTRACK.4AppOFZRA     <- rep(0, length(time.AppOFZRA))
  tempa.FLTRACK.4AppOFZRA   <- rep(0, length(time.AppOFZRA))
  tempd.FLTRACK.4AppOFZRA   <- rep(0, length(time.AppOFZRA))
  comnt.FLTRACK.4AppOFZRA   <- rep(0, length(time.AppOFZRA))
  for (i in 1:length(time.AppOFZRA)) {
    # find index of when chat time matches flight track time
    ind.FLTRACKTIME.4AppOFZRA[i] <- which(abs(difftime(flighttrack$ac.date.posix, time.AppOFZRA[i])) == min(abs(difftime(flighttrack$ac.date.posix, time.AppOFZRA[i])))) + 1
    # record time/lat/lon/alt/temp at AC location for time when IFI keyword was observed by flight scientist
    time.FLTRACK.4AppOFZRA[i]    <- flighttrack$ac.date.posix[ind.FLTRACKTIME.4AppOFZRA[i]]
    lat.FLTRACK.4AppOFZRA[i]     <- flighttrack$lat[ind.FLTRACKTIME.4AppOFZRA[i]]
    lon.FLTRACK.4AppOFZRA[i]     <- flighttrack$lon[ind.FLTRACKTIME.4AppOFZRA[i]]
    alt.FLTRACK.4AppOFZRA[i]     <- flighttrack$pres.alt.ft[ind.FLTRACKTIME.4AppOFZRA[i]]
    tempa.FLTRACK.4AppOFZRA[i]   <- flighttrack$amb.temp.degC[ind.FLTRACKTIME.4AppOFZRA[i]]
    tempd.FLTRACK.4AppOFZRA[i]   <- flighttrack$dew.temp.degC[ind.FLTRACKTIME.4AppOFZRA[i]]
    comnt.FLTRACK.4AppOFZRA[i]   <- as.character(chat.data$COMMENT[index.AppOFZRA.I2[i] - 1])
  }
}  else {
  time.FLTRACK.4AppOFZRA    <- -99
  lat.FLTRACK.4AppOFZRA     <- -99
  lon.FLTRACK.4AppOFZRA     <- -99
  alt.FLTRACK.4AppOFZRA     <- -99
  tempa.FLTRACK.4AppOFZRA   <- -99
  tempd.FLTRACK.4AppOFZRA   <- -99
  comnt.FLTRACK.4AppOFZRA   <- -99
}


# loop through all times when AppOFZDZ keyword was called out by flight scientist
if (length(time.AppOFZDZ) > 0) {
  ind.FLTRACKTIME.4AppOFZDZ <- rep(0, length(time.AppOFZDZ))
  time.FLTRACK.4AppOFZDZ    <- rep(0, length(time.AppOFZDZ))
  lat.FLTRACK.4AppOFZDZ     <- rep(0, length(time.AppOFZDZ))
  lon.FLTRACK.4AppOFZDZ     <- rep(0, length(time.AppOFZDZ))
  alt.FLTRACK.4AppOFZDZ     <- rep(0, length(time.AppOFZDZ))
  tempa.FLTRACK.4AppOFZDZ   <- rep(0, length(time.AppOFZDZ))
  tempd.FLTRACK.4AppOFZDZ   <- rep(0, length(time.AppOFZDZ))
  comnt.FLTRACK.4AppOFZDZ   <- rep(0, length(time.AppOFZDZ))
  for (i in 1:length(time.AppOFZDZ)) {
    # find index of when chat time matches flight track time
    ind.FLTRACKTIME.4AppOFZDZ[i] <- which(abs(difftime(flighttrack$ac.date.posix, time.AppOFZDZ[i])) == min(abs(difftime(flighttrack$ac.date.posix, time.AppOFZDZ[i])))) + 1
    # record time/lat/lon/alt/temp at AC location for time when IFI keyword was observed by flight scientist
    time.FLTRACK.4AppOFZDZ[i]    <- flighttrack$ac.date.posix[ind.FLTRACKTIME.4AppOFZDZ[i]]
    lat.FLTRACK.4AppOFZDZ[i]     <- flighttrack$lat[ind.FLTRACKTIME.4AppOFZDZ[i]]
    lon.FLTRACK.4AppOFZDZ[i]     <- flighttrack$lon[ind.FLTRACKTIME.4AppOFZDZ[i]]
    alt.FLTRACK.4AppOFZDZ[i]     <- flighttrack$pres.alt.ft[ind.FLTRACKTIME.4AppOFZDZ[i]]
    tempa.FLTRACK.4AppOFZDZ[i]   <- flighttrack$amb.temp.degC[ind.FLTRACKTIME.4AppOFZDZ[i]]
    tempd.FLTRACK.4AppOFZDZ[i]   <- flighttrack$dew.temp.degC[ind.FLTRACKTIME.4AppOFZDZ[i]]
    comnt.FLTRACK.4AppOFZDZ[i]   <- as.character(chat.data$COMMENT[index.AppOFZDZ.I2[i] - 1])
  }
}  else {
  time.FLTRACK.4AppOFZDZ    <- -99
  lat.FLTRACK.4AppOFZDZ     <- -99
  lon.FLTRACK.4AppOFZDZ     <- -99
  alt.FLTRACK.4AppOFZDZ     <- -99
  tempa.FLTRACK.4AppOFZDZ   <- -99
  tempd.FLTRACK.4AppOFZDZ   <- -99
  comnt.FLTRACK.4AppOFZDZ   <- -99
}

# loop through all times when ALLICE keyword was called out by flight scientist
if (length(time.ALLICE) > 0) {
  ind.FLTRACKTIME.4ALLICE <- rep(0, length(time.ALLICE))
  time.FLTRACK.4ALLICE    <- rep(0, length(time.ALLICE))
  lat.FLTRACK.4ALLICE     <- rep(0, length(time.ALLICE))
  lon.FLTRACK.4ALLICE     <- rep(0, length(time.ALLICE))
  alt.FLTRACK.4ALLICE     <- rep(0, length(time.ALLICE))
  tempa.FLTRACK.4ALLICE   <- rep(0, length(time.ALLICE))
  tempd.FLTRACK.4ALLICE   <- rep(0, length(time.ALLICE))
  comnt.FLTRACK.4ALLICE   <- rep(0, length(time.ALLICE))
  for (i in 1:length(time.ALLICE)) {
    # find index of when chat time matches flight track time
    ind.FLTRACKTIME.4ALLICE[i] <- which(abs(difftime(flighttrack$ac.date.posix, time.ALLICE[i])) == min(abs(difftime(flighttrack$ac.date.posix, time.ALLICE[i])))) + 1
    # record time/lat/lon/alt/temp at AC location for time when IFI keyword was observed by flight scientist
    time.FLTRACK.4ALLICE[i]    <- flighttrack$ac.date.posix[ind.FLTRACKTIME.4ALLICE[i]]
    lat.FLTRACK.4ALLICE[i]     <- flighttrack$lat[ind.FLTRACKTIME.4ALLICE[i]]
    lon.FLTRACK.4ALLICE[i]     <- flighttrack$lon[ind.FLTRACKTIME.4ALLICE[i]]
    alt.FLTRACK.4ALLICE[i]     <- flighttrack$pres.alt.ft[ind.FLTRACKTIME.4ALLICE[i]]
    tempa.FLTRACK.4ALLICE[i]   <- flighttrack$amb.temp.degC[ind.FLTRACKTIME.4ALLICE[i]]
    tempd.FLTRACK.4ALLICE[i]   <- flighttrack$dew.temp.degC[ind.FLTRACKTIME.4ALLICE[i]]
    comnt.FLTRACK.4ALLICE[i]   <- as.character(chat.data$COMMENT[index.ALLICE.I2[i] - 1])
  }
} else {
  time.FLTRACK.4ALLICE    <- -99
  lat.FLTRACK.4ALLICE     <- -99
  lon.FLTRACK.4ALLICE     <- -99
  alt.FLTRACK.4ALLICE     <- -99
  tempa.FLTRACK.4ALLICE   <- -99
  tempd.FLTRACK.4ALLICE   <- -99
  comnt.FLTRACK.4ALLICE   <- -99
}

# loop through all times when DENDS keyword was called out by flight scientist
if (length(time.DENDS) > 0) {
  ind.FLTRACKTIME.4DENDS <- rep(0, length(time.DENDS))
  time.FLTRACK.4DENDS    <- rep(0, length(time.DENDS))
  lat.FLTRACK.4DENDS     <- rep(0, length(time.DENDS))
  lon.FLTRACK.4DENDS     <- rep(0, length(time.DENDS))
  alt.FLTRACK.4DENDS     <- rep(0, length(time.DENDS))
  tempa.FLTRACK.4DENDS   <- rep(0, length(time.DENDS))
  tempd.FLTRACK.4DENDS   <- rep(0, length(time.DENDS))
  comnt.FLTRACK.4DENDS   <- rep(0, length(time.DENDS))
  for (i in 1:length(time.DENDS)) {
    # find index of when chat time matches flight track time
    ind.FLTRACKTIME.4DENDS[i] <- which(abs(difftime(flighttrack$ac.date.posix, time.DENDS[i])) == min(abs(difftime(flighttrack$ac.date.posix, time.DENDS[i])))) + 1
    # record time/lat/lon/alt/temp at AC location for time when IFI keyword was observed by flight scientist
    time.FLTRACK.4DENDS[i]    <- flighttrack$ac.date.posix[ind.FLTRACKTIME.4DENDS[i]]
    lat.FLTRACK.4DENDS[i]     <- flighttrack$lat[ind.FLTRACKTIME.4DENDS[i]]
    lon.FLTRACK.4DENDS[i]     <- flighttrack$lon[ind.FLTRACKTIME.4DENDS[i]]
    alt.FLTRACK.4DENDS[i]     <- flighttrack$pres.alt.ft[ind.FLTRACKTIME.4DENDS[i]]
    tempa.FLTRACK.4DENDS[i]   <- flighttrack$amb.temp.degC[ind.FLTRACKTIME.4DENDS[i]]
    tempd.FLTRACK.4DENDS[i]   <- flighttrack$dew.temp.degC[ind.FLTRACKTIME.4DENDS[i]]
    comnt.FLTRACK.4DENDS[i]   <- as.character(chat.data$COMMENT[index.DENDS.I2[i] - 1])
  }
}  else {
  time.FLTRACK.4DENDS    <- -99
  lat.FLTRACK.4DENDS     <- -99
  lon.FLTRACK.4DENDS     <- -99
  alt.FLTRACK.4DENDS     <- -99
  tempa.FLTRACK.4DENDS   <- -99
  tempd.FLTRACK.4DENDS   <- -99
  comnt.FLTRACK.4DENDS   <- -99
}

# loop through all times when PLATES keyword was called out by flight scientist
if (length(time.PLATES) > 0) {
  ind.FLTRACKTIME.4PLATES <- rep(0, length(time.PLATES))
  time.FLTRACK.4PLATES    <- rep(0, length(time.PLATES))
  lat.FLTRACK.4PLATES     <- rep(0, length(time.PLATES))
  lon.FLTRACK.4PLATES     <- rep(0, length(time.PLATES))
  alt.FLTRACK.4PLATES     <- rep(0, length(time.PLATES))
  tempa.FLTRACK.4PLATES   <- rep(0, length(time.PLATES))
  tempd.FLTRACK.4PLATES   <- rep(0, length(time.PLATES))
  comnt.FLTRACK.4PLATES   <- rep(0, length(time.PLATES))
  for (i in 1:length(time.PLATES)) {
    # find index of when chat time matches flight track time
    ind.FLTRACKTIME.4PLATES[i] <- which(abs(difftime(flighttrack$ac.date.posix, time.PLATES[i])) == min(abs(difftime(flighttrack$ac.date.posix, time.PLATES[i])))) + 1
    # record time/lat/lon/alt/temp at AC location for time when IFI keyword was observed by flight scientist
    time.FLTRACK.4PLATES[i]    <- flighttrack$ac.date.posix[ind.FLTRACKTIME.4PLATES[i]]
    lat.FLTRACK.4PLATES[i]     <- flighttrack$lat[ind.FLTRACKTIME.4PLATES[i]]
    lon.FLTRACK.4PLATES[i]     <- flighttrack$lon[ind.FLTRACKTIME.4PLATES[i]]
    alt.FLTRACK.4PLATES[i]     <- flighttrack$pres.alt.ft[ind.FLTRACKTIME.4PLATES[i]]
    tempa.FLTRACK.4PLATES[i]   <- flighttrack$amb.temp.degC[ind.FLTRACKTIME.4PLATES[i]]
    tempd.FLTRACK.4PLATES[i]   <- flighttrack$dew.temp.degC[ind.FLTRACKTIME.4PLATES[i]]
    comnt.FLTRACK.4PLATES[i]   <- as.character(chat.data$COMMENT[index.PLATES.I2[i] - 1])
  }
} else {
  time.FLTRACK.4PLATES    <- -99
  lat.FLTRACK.4PLATES     <- -99
  lon.FLTRACK.4PLATES     <- -99
  alt.FLTRACK.4PLATES     <- -99
  tempa.FLTRACK.4PLATES   <- -99
  tempd.FLTRACK.4PLATES   <- -99
  comnt.FLTRACK.4PLATES   <- -99
} 
#-------------------------------------------------

#-------------------------------------------------
# interim plotting of these results
#-------------------------------------------------
plot.new()
par(mfrow=c(3, 5))
plot(lon.FLTRACK.4MIXPHA,                               lat.FLTRACK.4MIXPHA,        main = "",             xlab = "LON",             ylab = "LAT",        pch = "*", lty=1, col="blue",                  panel.first = {grid()})
plot(lon.FLTRACK.4AppC,                                 lat.FLTRACK.4AppC,          main = "",             xlab = "LON",             ylab = "LAT",        pch = "*", lty=1, col="blue",                  panel.first = {grid()})
plot(lon.FLTRACK.4AppOFZRA,                             lat.FLTRACK.4AppOFZRA,      main = paste(as.character(icicle.2.df$date.posix[ind.fl.num]), " F", as.character(ICICLE.fl.num), sep=""), xlab = "LON",             ylab = "LAT",        pch = "*", lty=1, col="blue",                  panel.first = {grid()})
plot(lon.FLTRACK.4AppOFZDZ,                             lat.FLTRACK.4AppOFZDZ,      main = "",             xlab = "LON",             ylab = "LAT",        pch = "*", lty=1, col="blue",                  panel.first = {grid()})
plot(lon.FLTRACK.4ALLICE,                               lat.FLTRACK.4ALLICE,        main = "",             xlab = "LON",             ylab = "LAT",        pch = "*", lty=1, col="blue",                  panel.first = {grid()})
plot(tempa.FLTRACK.4MIXPHA,                             alt.FLTRACK.4MIXPHA/1000,   main = "MIXPHA",       xlab = "TEMP [degC]",     ylab = "ALT [kft]",  pch = "*", lty=1, col="blue", ylim = c(0, 20), xlim = c(-20, 5), panel.first = {grid()})
text(-5, 18, paste("N=", as.character(length(time.FLTRACK.4MIXPHA), sep="")))
plot(tempa.FLTRACK.4AppC,                               alt.FLTRACK.4AppC/1000,     main = "AppC",         xlab = "TEMP [degC]",     ylab = "ALT [kft]",  pch = "*", lty=1, col="blue", ylim = c(0, 20), xlim = c(-20, 5), panel.first = {grid()})
text(-5, 18, paste("N=", as.character(length(time.FLTRACK.4AppC), sep="")))
plot(tempa.FLTRACK.4AppOFZRA,                           alt.FLTRACK.4AppOFZRA/1000, main = "AppOFZRA",     xlab = "TEMP [degC]",     ylab = "ALT [kft]",  pch = "*", lty=1, col="blue", ylim = c(0, 20), xlim = c(-20, 5), panel.first = {grid()})
text(-5, 18, paste("N=", as.character(length(time.FLTRACK.4AppOFZRA), sep="")))
plot(tempa.FLTRACK.4AppOFZDZ,                           alt.FLTRACK.4AppOFZDZ/1000, main = "AppOFZDZ",     xlab = "TEMP [degC]",     ylab = "ALT [kft]",  pch = "*", lty=1, col="blue", ylim = c(0, 20), xlim = c(-20, 5), panel.first = {grid()})
text(-5, 18, paste("N=", as.character(length(time.FLTRACK.4AppOFZDZ), sep="")))
plot(tempa.FLTRACK.4ALLICE,                             alt.FLTRACK.4ALLICE/1000,   main = "ALLICE",       xlab = "TEMP [degC]",     ylab = "ALT [kft]",  pch = "*", lty=1, col="blue", ylim = c(0, 20), xlim = c(-20, 5), panel.first = {grid()})
text(-5, 18, paste("N=", as.character(length(time.FLTRACK.4ALLICE), sep="")))
plot((tempa.FLTRACK.4MIXPHA-tempd.FLTRACK.4MIXPHA),     alt.FLTRACK.4MIXPHA/1000,   main = "",             xlab = "DEW-DEPR [degC]", ylab = "ALT [kft]",  pch = "*", lty=1, col="blue", ylim = c(0, 20), xlim = c(-5, 5), panel.first = {grid()})
plot((tempa.FLTRACK.4AppC-tempd.FLTRACK.4AppC),         alt.FLTRACK.4AppC/1000,     main = "",             xlab = "DEW-DEPR [degC]", ylab = "ALT [kft]",  pch = "*", lty=1, col="blue", ylim = c(0, 20), xlim = c(-5, 5), panel.first = {grid()})
plot((tempa.FLTRACK.4AppOFZRA-tempd.FLTRACK.4AppOFZRA), alt.FLTRACK.4AppOFZRA/1000, main = "",             xlab = "DEW-DEPR [degC]", ylab = "ALT [kft]",  pch = "*", lty=1, col="blue", ylim = c(0, 20), xlim = c(-5, 5), panel.first = {grid()})
plot((tempa.FLTRACK.4AppOFZDZ-tempd.FLTRACK.4AppOFZDZ), alt.FLTRACK.4AppOFZDZ/1000, main = "",             xlab = "DEW-DEPR [degC]", ylab = "ALT [kft]",  pch = "*", lty=1, col="blue", ylim = c(0, 20), xlim = c(-5, 5), panel.first = {grid()})
plot((tempa.FLTRACK.4ALLICE-tempd.FLTRACK.4ALLICE),     alt.FLTRACK.4ALLICE/1000,   main = "",             xlab = "DEW-DEPR [degC]", ylab = "ALT [kft]",  pch = "*", lty=1, col="blue", ylim = c(0, 20), xlim = c(-5, 5), panel.first = {grid()})

#-------------------------------------------------
# save off image
#-------------------------------------------------
if (output.plot.flag == 1) {
  dev.copy(png, paste(output.plot.dir, as.character(icicle.2.df$date.posix[ind.fl.num]), "_F", as.character(ICICLE.fl.num), "_autosort_keyword.png", sep = ""))
  dev.off()
}
#-------------------------------------------------

#-------------------------------------------------
# save off flight track time, alt, lat, lon, heading information
#-------------------------------------------------
# first create data frames with the values for each category
if (output.data.flag == 1) {
  
  if (length(time.MIXPHA) > 0) {
    MIXPHA.df   <- data.frame(time.MIXPHA,   lat.FLTRACK.4MIXPHA,   lon.FLTRACK.4MIXPHA,   alt.FLTRACK.4MIXPHA,   tempa.FLTRACK.4MIXPHA,   tempd.FLTRACK.4MIXPHA, comnt.FLTRACK.4MIXPHA)
    write.csv(MIXPHA.df,   file = paste(output.data.dir, as.character(icicle.2.df$date.posix[ind.fl.num]), "_F", as.character(ICICLE.fl.num), "_MIXPHA_autosort_keyword.csv",   sep=""))
  }
  if (length(time.AppC) > 0) {
    AppC.df     <- data.frame(time.AppC,     lat.FLTRACK.4AppC,     lon.FLTRACK.4AppC,     alt.FLTRACK.4AppC,     tempa.FLTRACK.4AppC,     tempd.FLTRACK.4AppC, comnt.FLTRACK.4AppC)
    write.csv(AppC.df,     file = paste(output.data.dir, as.character(icicle.2.df$date.posix[ind.fl.num]), "_F", as.character(ICICLE.fl.num), "_AppC_autosort_keyword.csv",     sep=""))
  }
  if (length(time.AppOFZRA) > 0) {
    AppOFZRA.df <- data.frame(time.AppOFZRA, lat.FLTRACK.4AppOFZRA, lon.FLTRACK.4AppOFZRA, alt.FLTRACK.4AppOFZRA, tempa.FLTRACK.4AppOFZRA, tempd.FLTRACK.4AppOFZRA, comnt.FLTRACK.4AppOFZRA)
    write.csv(AppOFZRA.df, file = paste(output.data.dir, as.character(icicle.2.df$date.posix[ind.fl.num]), "_F", as.character(ICICLE.fl.num), "_AppOFZRA_autosort_keyword.csv", sep=""))
  }
  if (length(time.AppOFZDZ) > 0) {
    AppOFZDZ.df <- data.frame(time.AppOFZDZ, lat.FLTRACK.4AppOFZDZ, lon.FLTRACK.4AppOFZDZ, alt.FLTRACK.4AppOFZDZ, tempa.FLTRACK.4AppOFZDZ, tempd.FLTRACK.4AppOFZDZ, comnt.FLTRACK.4AppOFZDZ)
    write.csv(AppOFZDZ.df, file = paste(output.data.dir, as.character(icicle.2.df$date.posix[ind.fl.num]), "_F", as.character(ICICLE.fl.num), "_AppOFZDZ_autosort_keyword.csv", sep=""))
  }
  if (length(time.ALLICE) > 0) {
    ALLICE.df   <- data.frame(time.ALLICE,   lat.FLTRACK.4ALLICE,   lon.FLTRACK.4ALLICE,   alt.FLTRACK.4ALLICE,   tempa.FLTRACK.4ALLICE,   tempd.FLTRACK.4ALLICE, comnt.FLTRACK.4ALLICE)
    write.csv(ALLICE.df,   file = paste(output.data.dir, as.character(icicle.2.df$date.posix[ind.fl.num]), "_F", as.character(ICICLE.fl.num), "_ALLICE_autosort_keyword.csv",   sep=""))
  }
  if (length(time.DENDS) > 0) {
    DENDS.df    <- data.frame(time.DENDS,   lat.FLTRACK.4DENDS,   lon.FLTRACK.4DENDS,   alt.FLTRACK.4DENDS,   tempa.FLTRACK.4DENDS,   tempd.FLTRACK.4DENDS, comnt.FLTRACK.4DENDS)
    write.csv(DENDS.df,   file = paste(output.data.dir, as.character(icicle.2.df$date.posix[ind.fl.num]), "_F", as.character(ICICLE.fl.num), "_DENDS_autosort_keyword.csv",   sep=""))
  }
  if (length(time.PLATES) > 0) {
    PLATES.df   <- data.frame(time.PLATES,   lat.FLTRACK.4PLATES, lon.FLTRACK.4PLATES,   alt.FLTRACK.4PLATES,   tempa.FLTRACK.4PLATES,   tempd.FLTRACK.4PLATES, comnt.FLTRACK.4PLATES)
    write.csv(PLATES.df,   file = paste(output.data.dir, as.character(icicle.2.df$date.posix[ind.fl.num]), "_F", as.character(ICICLE.fl.num), "_PLATES_autosort_keyword.csv",   sep=""))
  }
} # end of if (output.data.flag == 1)
#-------------------------------------------------
# load appropriate RadIA interest volume netcdf file
#-------------------------------------------------

#-------------------------------------------------
# find the radar tilt that results in closest altitude match to range/azim from radar that ac was at during indexed time
#-------------------------------------------------

#-------------------------------------------------
# locate/save off radia ints for small/FZDZ/MIXPHA 
#-------------------------------------------------

# repeat this for each ICICLE Flight where RadIA data had been processed in realtime