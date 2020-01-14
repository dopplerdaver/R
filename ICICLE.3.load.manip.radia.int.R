#---------------------------------
# Name:     ICICLE.3.load.manip.radia.int.R
#
# Purpose:  NCAR scientists won't have even quicklook data until at least september 2019.
#           In order to make progress with algorithm skill analysis, the following process is conducted:
#
#           Execute the following scripts in the following order:
#
#           1. ICICLE.1.related.data.frames.R
#               a. builds 'NEXRAD.site.df', 'icicle.df' and 'icicle.chat.ificat.df' data frames
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
# Created:  3.25.2019 dserke
#---------------------------------

#---------------------------------
# user input
#---------------------------------
#   user.yyyymmdd must exist
user.yyyymmdd                     <- "20190227" # ICICLE F24
user.hhmm                         <- "0055"
user.hhmm.numeric                 <- as.numeric(user.hhmm)

#---------------------------------
# define directories and filelists in those dirs
#---------------------------------

SRPC.mom.dir                <-
SRCC.mom.dir                <-
MOCC.mom.dir                <-

SRPC.RadIA.INT.dir          <- file.path(paste("/d1/serke/projects/case_studies/ICICLE_2018/data/nc_radia/KDVN/polar/", user.yyyymmdd, "/", sep = "")) # ICICLE Flights
SRCC.RadIA.INT.dir          <- file.path(paste("/d1/serke/projects/case_studies/ICICLE_2018/data/nc_radia/KDVN/cart/",  user.yyyymmdd, "/", sep = "")) # ICICLE Flights
MOCC.RadIA.INT.dir          <- file.path(paste("/d1/serke/projects/case_studies/ICICLE_2018/data/nc_radia/KDVN/cart/",  user.yyyymmdd, "/", sep = "")) # ICICLE Flights

# for radia, naming convention: SRPC=singrad,polcoord, SRCC=singrad,cartcoord, MOCC=mosaic,cartcoord
SRPC.RadIA.INT.filelist     <- dir(SRPC.RadIA.INT.dir,  pattern = "[.nc]",       full.names = FALSE, ignore.case = TRUE)
SRPC.RadIA.INT.gz.filelist  <- dir(SRPC.RadIA.INT.dir,  pattern = "[.nc.gz]",    full.names = FALSE, ignore.case = TRUE)

SRCC.RadIA.INT.filelist     <- dir(SRCC.RadIA.INT.dir,  pattern = "[.nc]",       full.names = FALSE, ignore.case = TRUE)
MOCC.RadIA.INT.filelist     <- dir(MOCC.RadIA.INT.dir,  pattern = "[.nc]",       full.names = FALSE, ignore.case = TRUE)

#-------------------------------------------------
# load appropriate RadIA interest volume netcdf file
#-------------------------------------------------
#   find the RadIA nc file closest to user specified hhmm
ind.min.time.diff           <- which.min(abs(as.numeric(substr(SRPC.RadIA.INT.gz.filelist, 1, 4)) - user.hhmm.numeric))

#   unzip nc file if its currently gzipped
system(paste('gunzip ', SRPC.RadIA.INT.dir, SRPC.RadIA.INT.gz.filelist[ind.min.time.diff], sep =""))

#

#-------------------------------------------------
# find the radar tilt that results in closest altitude match to range/azim from radar that ac was at during indexed time
#-------------------------------------------------

#-------------------------------------------------
# locate/save off radia ints for small/FZDZ/MIXPHA 
#-------------------------------------------------

# repeat this for each ICICLE Flight where RadIA data had been processed in realtime