
#-------------------------------------------------
# Script Name:  ICICLE.1.related.data.frames.R
#
# Purpose:      builds 'NEXRAD.site.df', icicle.df' and 'icicle.chat.ificat.df' data frames
#
#            	Execute the following scripts in the following order:
#
#           	1. ICICLE.1.related.data.frames.R
#               	a. builds 'NEXRAD.site.df', icicle.df' and 'icicle.chat.ificat.df' data frames
#           	2. ICICLE.2.sort.chat.log.R
#           		a. Search through CNV flight chat log to find specified keyword occurrences, then find date/time of occurrences
#           		b. If keyword occurrence is within primary radar domain (within ~100km), search CNV flight track log
#              	 		 for the lat/lon/alt/tempa/tempd values.  
#               	c. Plot these metadata for individual flights and whole campaign
#               	d. Save off csv files of these data from every flight
#          	 3. ICICLE.3.load.manip.radia.int.R
#              		a. Find closest radar volume in time and find the closest RadIA INTs to the keyword occurrences.
#              		b. Plot statistics on num occur, alts, mean/median/boxwhisker plots of INT values, etc
#
#
#           	Approx number of hours of data to assess: 27hr and 43min plus primary time from F7,17,18,23, & 26
#           	Number of flights to assess:              26 research flights, 19 with processed RadIA
#
# Created:      1.30.2019 dserke
#-------------------------------------------------

library(gtable)
library(grid)
library(gridExtra)

# set up dir path to NEXRAD site data file
nexrad.site.dataframetxt.dir    <- file.path("/d1/serke/projects/case_studies/SNOWIE/data/RadIA_data/nexrad_site_data/")
# read in the NEXRAD site location text file
NEXRAD.site.df                  <- read.csv(paste(nexrad.site.dataframetxt.dir, "nexrad_site.csv", sep = ""), header = FALSE, sep = ",", dec = ".", stringsAsFactors=FALSE)
colnames(NEXRAD.site.df)        <- c("NCDCID", "ICAO", "WBAN", "radname", "COUNTRY", "STATE", "COUNTY", "lat", "lon", "elev", "GMTdiff", "STN_TYPE")
#head(NEXRAD.site.df)

# constants
fl.hours.avail                  <- 120 # alloted flight hours from FAA
camp.tot.days                   <- 40  # num days between 1-27 and 3-9 during 2019

numflights.total                <- 29
numflights.research             <- 26
numflights.woRadIA              <- 7
numflights.wRadIA               <- 19

MIXPHA.numflights.wRadIA        <- 16
AppC.numflights.wRadIA          <- 7
AppOFZRA.numflights.wRadIA      <- 4
AppOFZDZ.numflights.wRadIA      <- 14
ALLICE.numflights.wRadIA        <- 12
DENDS.numflights.wRadIA         <- 6  
PLATES.numflights.wRadIA        <- 2  

MIXPHA.numchatrefs.wRadIA       <- 62
AppC.numchatrefs.wRadIA         <- 20
AppOFZRA.numchatrefs.wRadIA     <- 5
AppOFZDZ.numchatrefs.wRadIA     <- 162
ALLICE.numchatrefs.wRadIA       <- 30
DENDS.numchatrefs.wRadIA        <- 18 
PLATES.numchatrefs.wRadIA       <- 2  

# build fields for icicle.chat.ificat.df
numFLTS                         <- c(MIXPHA.numflights.wRadIA,  AppC.numflights.wRadIA,  AppOFZRA.numflights.wRadIA,  AppOFZDZ.numflights.wRadIA,  ALLICE.numflights.wRadIA,  DENDS.numflights.wRadIA,  PLATES.numflights.wRadIA, numflights.wRadIA, numflights.research, numflights.total)
numCHAT                         <- c(MIXPHA.numchatrefs.wRadIA, AppC.numchatrefs.wRadIA, AppOFZRA.numchatrefs.wRadIA, AppOFZDZ.numchatrefs.wRadIA, ALLICE.numchatrefs.wRadIA, DENDS.numchatrefs.wRadIA, PLATES.numchatrefs.wRadIA, NA, NA, NA)
icicle.chat.ificat.df           <- data.frame(numFLTS, numCHAT)
colnames(icicle.chat.ificat.df) <- c("Num-Flights", "Num-chat-refs")
rownames(icicle.chat.ificat.df) <- c("MIX-PHA", "App-C", "App-O FZRA", "App-O FZDZ", "ALL-ICE", "DENDRITES", "PLATES", "W/PROC-RadIA", "W/RESEARCH", "TOTAL")
head(icicle.chat.ificat.df)

# build fields
date.posix            <- c("2019-01-27", "2019-01-28", "2019-01-29", "2019-01-30", "2019-01-31", "2019-02-01", "2019-02-02", "2019-02-03", "2019-02-04", "2019-02-05", "2019-02-06", "2019-02-07", "2019-02-08", "2019-02-09", "2019-02-10", "2019-02-11", "2019-02-12", "2019-02-12", "2019-02-13", "2019-02-14", "2019-02-15", "2019-02-16", "2019-02-17", "2019-02-17", "2019-02-18", "2019-02-19", "2019-02-20", "2019-02-21", "2019-02-22", "2019-02-23", "2019-02-24", "2019-02-24", "2019-02-25", "2019-02-26", "2019-02-26", "2019-02-27", "2019-02-28", "2019-03-01", "2019-03-02", "2019-03-02", "2019-03-03", "2019-03-04", "2019-03-05", "2019-03-06", "2019-03-07") # date.posix
camp.day.num          <- c(           1,            2,            3,            4,            5,            6,            7,            8,            9,           10,           11,           12,           13,           14,           15,           16,           17,           17,           18,           19,           20,           21,           21,           22,           23,           24,           25,           26,           27,           28,           29,           29,           30,           31,           31,           32,           33,           34,           35,           35,           36,           37,           38,           39,           40) # camp.day.num
case.hours.used       <- c(         3.5,          3.9,          4.0,            0,          4.0,            0,            0,            0,          4.0,          4.3,          4.0,          4.0,            3,            0,            0,          3.5,          4.0,          3.0,            0,          3.8,          4.3,          4.5,          4.2,          4.0,            0,            0,            0,            0,            5,          4.0,          3.9,          4.0,            0,          3.0,          4.0,            0,          3.5,           NA,          4.3,          3.0,            0,            0,          4.5,            0,           NA) # case.hrs.used
fl.num                <- c(           2,            3,            4,           NA,            5,           NA,           NA,           NA,            6,            7,            8,            9,           10,           NA,           NA,           11,           12,           13,           NA,           14,           15,           16,           17,           18,           NA,           NA,           NA,           NA,           19,           20,           21,           22,           NA,           23,           24,           NA,           25,           NA,           26,           27,           NA,           NA,           28,           NA,           29) # fl.num
fl.srt                <- c(     "18:00",      "20:05",      "18:00",           NA,      "16:00",           NA,           NA,           NA,      "12:23",      "17:05",      "19:02",      "18:00",      "14:00",           NA,           NA,      "15:00",      "11:07",      "16:55",           NA,      "16:19",      "19:34",      "01:16",      "12:04",      "16:55",           NA,           NA,           NA,           NA,      "16:24",      "12:03",      "11:56",      "17:07",           NA,      "17:24",      "22:42",           NA,      "17:55",           NA,      "10:32",      "16:11",           NA,           NA,      "11:54",           NA,      "16:16") # fl.srt
fl.end                <- c(     "21:30",      "00:00",      "22:00",           NA,      "20:00",           NA,           NA,           NA,      "16:23",      "22:00",      "22:52",      "22:00",      "17:00",           NA,           NA,      "18:25",      "15:08",      "20:00",           NA,      "20:37",      "23:59",      "05:44",      "16:06",      "20:17",           NA,           NA,           NA,           NA,      "21:30",      "16:07",      "15:49",      "21:06",           NA,      "20:51",      "02:44",           NA,      "21:20",           NA,      "14:51",      "19:11",           NA,           NA,      "16:27",           NA,      "21:06") # fl.end
radia.priority.1to5   <- c(          NA,            2,            4,           NA,            3,           NA,           NA,           NA,            2,            2,            1,            1,            5,           NA,           NA,            4,            2,            1,           NA,            2,            2,            1,            2,            4,           NA,           NA,           NA,           NA,            5,            2,            1,            2,           NA,            5,            1,           NA,            2,           NA,            1,            3,           NA,           NA,            2,           NA,            2) #
IFI.cond              <- c(  "FERRYOTT",  "C/SLD/FDZ","C/CV/LES/cold/MP",      NA,   "C/MP/ICE",           NA,           NA,           NA, "ZBB/FRA/MP", "CV/FRA/SEV", "FDZ/SLD/UP", "MP/FDZ/FRA",   "FERRYOTT",           NA,           NA,   "FERRYOTT","C/FRA/FDZ/MP",  "FDZ2C/MP",           NA,   "C/MP/ICE",  "FRA/MP/IP", "FRA/FDZ/MP", "FDZ/SLD/MP",      "C/FDZ",           NA,           NA,           NA,           NA,    "CAL/FDZ","FDZ/SLD/FRA/HI",  "FDZ/VC",   "C/ICE2MP",           NA,       "C/CA","SLDBEST/FDZ",           NA,  "ALL/trans",           NA,"FZDZ/MIXPHA","endFZ/LE-SCu",          NA,           NA,"LES/cold/C/SCu/MP",     NA,        "C/0") #
pri.radar.name        <- c(          NA,       "KGRR",       "KGRR",           NA,       "KDMX",           NA,           NA,           NA,       "KMQT",       "KDVN",       "KDMX",       "KGRB",           NA,           NA,           NA,       "KMKX",       "KGRR",       "KGRR",           NA,       "KDVN",       "KPAH",       "KPAH",       "KDVN",       "KILX",           NA,           NA,           NA,           NA,       "KDVN",       "KMKX",       "KIWX",       "KGRR",           NA,       "KDVN",       "KDVN",           NA,       "KVWX",           NA,       "KGRR",       "KIWX",           NA,           NA,       "KGRR",           NA,       "KDMX") # prime.radar.name
pri.radar.srt         <- c(          NA,      "20:30",      "18:37",           NA,      "18:10",           NA,           NA,           NA,      "13:21",      "20:00",      "20:00",      "18:40",           NA,           NA,           NA,      "18:12",      "11:42",      "17:32",           NA,      "16:31",      "20:37",      "01:54",      "12:13",      "17:41",           NA,           NA,           NA,           NA,      "18:05",      "12:10",      "12:32",      "17:38",           NA,      "17:46",      "23:01",           NA,      "19:00",           NA,      "11:15",      "16:11",           NA,           NA,      "12:30",           NA,      "17:01") #
pri.radar.end         <- c(          NA,      "23:00",      "19:20",           NA,      "19:28",           NA,           NA,           NA,      "15:10",      "21:00",      "22:13",      "21:05",           NA,           NA,           NA,      "18:25",      "14:07",      "18:35",           NA,      "18:16",      "23:15",      "03:18",      "13:53",      "18:48",           NA,           NA,           NA,           NA,      "18:40",      "15:00",      "14:32",      "20:16",           NA,      "18:06",      "02:25",           NA,      "21:05",           NA,      "12:09",      "17:04",           NA,           NA,      "15:17",           NA,      "19:51") #
pri.radar.pcd         <- c(          NA,            0,            0,           NA,            0,           NA,           NA,           NA,            0,          1.0,          0.0,          1.0,           NA,           NA,           NA,          1.0,          1.0,          1.0,           NA,          1.0,            0,            0,          1.0,          1.0,           NA,           NA,           NA,           NA,            1,            1,            1,            1,           NA,            1,            1,           NA,            1,           NA,            1,            1,           NA,           NA,            1,           NA,            1) #
sec.radar.name        <- c(          NA,       "KLOT",       "KMKW",           NA,       "KDVN",           NA,           NA,           NA,       "KGRB",       "KARX",       "KDVN",       "KMKX",           NA,           NA,           NA,       "KGRR",       "KMKX",       "KMKX",           NA,       "KARX",       "KVWX",       "KSGF",       "KILX",       "KLSX",           NA,           NA,           NA,           NA,         "NA",       "KARX",       "KGRR",       "KIWX",           NA,       "KLOT",       "KDMX",           NA,       "KILX",           NA,       "KIWX",       "KGRR",           NA,           NA,         "NA",           NA,       "KLSX") #
sec.radar.srt         <- c(          NA,      "20:00",      "19:25",           NA,      "16:56",           NA,           NA,           NA,      "15:18",           NA,           NA,      "18:10",           NA,           NA,           NA,      "17:20",           NA,      "17:05",           NA,      "17:25",      "23:18",      "03:46",      "14:07",      "18:55",           NA,           NA,           NA,           NA,           NA,      "15:12",      "14:45",      "17:07",           NA,           NA,      "22:42",           NA,      "21:22",           NA,      "12:16",      "17:40",           NA,           NA,           NA,           NA,           NA) #
sec.radar.end         <- c(          NA,      "20:30",      "20:00",           NA,      "18:07",           NA,           NA,           NA,      "15:40",           NA,           NA,      "18:46",           NA,           NA,           NA,      "17:50",           NA,      "17:20",           NA,      "19:15",      "23:59",      "03:56",      "15:40",      "19:16",           NA,           NA,           NA,           NA,           NA,      "15:12",      "15:15",      "17:30",           NA,           NA,      "02:44",           NA,      "21:45",           NA,      "14:51",           NA,           NA,           NA,           NA,           NA,           NA) #
sec.radar.pcd         <- c(          NA,            0,            0,           NA,            0,           NA,           NA,           NA,            0,            1,            0,            1,           NA,           NA,           NA,            0,            1,            1,           NA,            1,            0,            0,            0,            1,           NA,           NA,           NA,           NA,           NA,            1,            1,            1,           NA,            1,            1,           NA,            1,           NA,            1,            1,           NA,           NA,            1,           NA,            1) #
radia.rpt.done        <- c(          NA,            1,            1,           NA,            1,           NA,           NA,           NA,          0.8,          0.8,          0.8,            1,           NA,           NA,           NA,          0.9,          0.9,          0.9,           NA,          0.9,          0.7,          0.7,          1.0,          1.0,           NA,           NA,           NA,           NA,          0.9,          1.0,          0.7,          0.8,           NA,          0.9,          0.9,           NA,          1.0,           NA,          0.8,          1.0,           NA,           NA,          1.0,           NA,          0.5) #
wx.rpt.resp           <- c(          NA,           NA,           NA,           NA,           NA,           NA,           NA,           NA,           NA,           NA,           NA,           NA,           NA,           NA,           NA,           NA,           NA,           NA,           NA,           NA,           NA,           NA,           NA,           NA,           NA,           NA,           NA,           NA,           NA,           NA,           NA,           NA,           NA,           NA,           NA,            1,            1,            1,            1,           NA,            1,            1,            1,            1,            1) #
wx.rpt.done           <- c(          NA,           NA,           NA,           NA,           NA,           NA,           NA,           NA,           NA,           NA,           NA,           NA,           NA,           NA,           NA,           NA,           NA,           NA,           NA,           NA,           NA,           NA,           NA,           NA,           NA,           NA,           NA,           NA,           NA,           NA,           NA,           NA,           NA,           NA,           NA,            1,            0,            1,            1,           NA,            1,            1,            1,            1,            1) #

tot.hours.used        <- cumsum(case.hours.used) 
tot.hours.ave         <- rep(0, length(camp.day.num))
for (i in 1 : length(camp.day.num)) {
  tot.hours.ave[i] <- (fl.hours.avail / camp.tot.days) * i
}
tot.hours.ave.new <- c(as.numeric(tot.hours.ave[1:17]), as.numeric(tot.hours.ave[17]), as.numeric(tot.hours.ave[18:19]), as.numeric(tot.hours.ave[20]), as.numeric(tot.hours.ave[20:29]), as.numeric(tot.hours.ave[29]), as.numeric(tot.hours.ave[30:33]), as.numeric(tot.hours.ave[33]), as.numeric(tot.hours.ave[34:length(tot.hours.ave)]))

# build data frame
icicle.df             <- data.frame(date.posix, camp.day.num, fl.num, fl.srt, fl.end, case.hours.used, tot.hours.used, pri.radar.name, pri.radar.srt, pri.radar.end, pri.radar.pcd, sec.radar.name, sec.radar.srt, sec.radar.end, sec.radar.pcd, radia.priority.1to5, radia.rpt.done, wx.rpt.resp, wx.rpt.done)
icicle.table.df       <- data.frame(date.posix, fl.num, IFI.cond, radia.priority.1to5, pri.radar.name, pri.radar.pcd, sec.radar.name, sec.radar.pcd, radia.rpt.done)
icicle.2.df           <- data.frame(date.posix, fl.num, IFI.cond, radia.priority.1to5, pri.radar.name, pri.radar.pcd, pri.radar.srt, pri.radar.end, sec.radar.name, sec.radar.pcd)
# info about data frame
head(icicle.df)
str(icicle.df)

# plotting

# plot timeseries of flight hours used (black w/ dots) versus mean hours used (red line)
plot.new()
par(mfrow = c(1, 1))
#by posix date
#plot(icicle.df$date.posix, icicle.df$tot.hours.ave, type = "p", col = "red", xlab = "day num", ylab = "tot hrs used", main = "ICICLE hrs used")
#by day num sequential
plot(icicle.df$camp.day.num, icicle.df$tot.hours.ave, type = "l", col = "red", xlab = "day num", ylab = "tot hrs used", main = "ICICLE hrs used")
lines(icicle.df$camp.day.num, icicle.df$tot.hours.used, type = "o", col = "black")
#lines(lowess(icicle.df$tot.hours.used), col="blue")
# TO DO: PLOT fl.num ABOVE EACH BLACK CIRCLE ON PLOT
legend(30, 60, c("hours ave", "hours used"), col = c(1, 2), cex = 1.0)
grid()

# plot table of number of flights/keyword mentions in chat for each IFI category 
fieldsfortable.1 <- icicle.chat.ificat.df
gridfortable.1   <- tableGrob(fieldsfortable.1, theme = ttheme_default(base_size=14, base_colour="black", base_family="", parse=FALSE, padding=unit(c(2, 2), "mm")))
gridfortable.1   <- gtable_add_grob(gridfortable.1, grobs = rectGrob(gp = gpar(fill = NA, lwd = 2)), t = 2, b = nrow(gridfortable.1), l = 1, r = ncol(gridfortable.1))
gridfortable.1   <- gtable_add_grob(gridfortable.1, grobs = rectGrob(gp = gpar(fill = NA, lwd = 2)), t = 1, l = 1, r = ncol(gridfortable.1))
gridfortable.1   <- gtable_add_grob(gridfortable.1, grobs = segmentsGrob( x0 = unit(0,"npc"), y0 = unit(0,"npc"), x1 = unit(1,"npc"), y1 = unit(0,"npc"), gp = gpar(lwd = 2.0)), t = 8, b = 3, l = 3, r = 3)
gridfortable.1   <- gtable_add_grob(gridfortable.1, grobs = segmentsGrob( x0 = unit(0,"npc"), y0 = unit(0,"npc"), x1 = unit(1,"npc"), y1 = unit(0,"npc"), gp = gpar(lwd = 2.0)), t = 8, b = 3, l = 2, r = 3)
gridfortable.1   <- gtable_add_grob(gridfortable.1, grobs = segmentsGrob( x0 = unit(0,"npc"), y0 = unit(0,"npc"), x1 = unit(1,"npc"), y1 = unit(0,"npc"), gp = gpar(lwd = 2.0)), t = 8, b = 3, l = 1, r = 3)
grid.newpage()
grid.draw(gridfortable.1)

# plot table of which flights have been/still require RadIA processing
fieldsfortable.2 <- icicle.table.df
gridfortable.2   <- tableGrob(fieldsfortable.2, theme = ttheme_default(base_size=8, base_colour="black", base_family="", parse=FALSE, padding=unit(c(0.5, 2), "mm")))
gridfortable.2   <- gtable_add_grob(gridfortable.2, grobs = rectGrob(gp = gpar(fill = NA, lwd = 2)), t = 2, b = nrow(gridfortable.2), l = 1, r = ncol(gridfortable.2))
gridfortable.2   <- gtable_add_grob(gridfortable.2, grobs = rectGrob(gp = gpar(fill = NA, lwd = 2)), t = 1, l = 1, r = ncol(gridfortable.2))
find_cell <- function(table, row, col, name="core-fg"){
  l <- table$layout
  which(l$t==row & l$l==col & l$name==name)
}
# primary radar unprocessed indices
ind01                                  <- find_cell(gridfortable.2,  3, 7, "core-bg")
ind02                                  <- find_cell(gridfortable.2,  4, 7, "core-bg")
ind03                                  <- find_cell(gridfortable.2,  6, 7, "core-bg")
ind04                                  <- find_cell(gridfortable.2, 10, 7, "core-bg")
ind05                                  <- find_cell(gridfortable.2, 12, 7, "core-bg")
ind06                                  <- find_cell(gridfortable.2, 22, 7, "core-bg")
ind07                                  <- find_cell(gridfortable.2, 23, 7, "core-bg")
# top radar priority cases indices
ind21                                  <- find_cell(gridfortable.2,  3, 9, "core-bg")
ind22                                  <- find_cell(gridfortable.2,  4, 9, "core-bg")
ind23                                  <- find_cell(gridfortable.2,  6, 9, "core-bg")
ind24                                  <- find_cell(gridfortable.2, 10, 9, "core-bg")
ind25                                  <- find_cell(gridfortable.2, 12, 9, "core-bg")
ind26                                  <- find_cell(gridfortable.2, 22, 9, "core-bg")
ind27                                  <- find_cell(gridfortable.2, 23, 9, "core-bg")
ind28                                  <- find_cell(gridfortable.2, 24, 9, "core-bg")
# secondary radar unprocessed indices
ind11                                  <- find_cell(gridfortable.2, 13, 5, "core-bg")
ind12                                  <- find_cell(gridfortable.2, 12, 5, "core-bg")
ind13                                  <- find_cell(gridfortable.2, 19, 5, "core-bg")
ind14                                  <- find_cell(gridfortable.2, 23, 5, "core-bg")
ind15                                  <- find_cell(gridfortable.2, 32, 5, "core-bg")
ind16                                  <- find_cell(gridfortable.2, 36, 5, "core-bg")
ind17                                  <- find_cell(gridfortable.2, 40, 5, "core-bg")
# FZRA cases indices
ind31                                  <- find_cell(gridfortable.2, 10, 6, "core-bg")
ind32                                  <- find_cell(gridfortable.2, 11, 6, "core-bg")
ind33                                  <- find_cell(gridfortable.2, 13, 6, "core-bg")
ind34                                  <- find_cell(gridfortable.2, 18, 6, "core-bg")
ind35                                  <- find_cell(gridfortable.2, 19, 6, "core-bg")
ind36                                  <- find_cell(gridfortable.2, 22, 6, "core-bg")
ind37                                  <- find_cell(gridfortable.2, 23, 6, "core-bg")
ind38                                  <- find_cell(gridfortable.2, 31, 6, "core-bg")
ind39                                  <- find_cell(gridfortable.2, 38, 6, "core-bg")
#
gridfortable.2$grobs[ind01][[1]][["gp"]] <- gpar(fontsize = 15)
gridfortable.2$grobs[ind02][[1]][["gp"]] <- gpar(fontsize = 15)
gridfortable.2$grobs[ind03][[1]][["gp"]] <- gpar(fontsize = 15)
gridfortable.2$grobs[ind04][[1]][["gp"]] <- gpar(fontsize = 15)
gridfortable.2$grobs[ind05][[1]][["gp"]] <- gpar(fontsize = 15)
gridfortable.2$grobs[ind06][[1]][["gp"]] <- gpar(fontsize = 15)
gridfortable.2$grobs[ind07][[1]][["gp"]] <- gpar(fontsize = 15)
gridfortable.2$grobs[ind21][[1]][["gp"]] <- gpar(fontsize = 15)
gridfortable.2$grobs[ind22][[1]][["gp"]] <- gpar(fontsize = 15)
gridfortable.2$grobs[ind23][[1]][["gp"]] <- gpar(fontsize = 15)
gridfortable.2$grobs[ind24][[1]][["gp"]] <- gpar(fontsize = 15)
gridfortable.2$grobs[ind25][[1]][["gp"]] <- gpar(fontsize = 15)
gridfortable.2$grobs[ind26][[1]][["gp"]] <- gpar(fontsize = 15)
gridfortable.2$grobs[ind27][[1]][["gp"]] <- gpar(fontsize = 15)
gridfortable.2$grobs[ind28][[1]][["gp"]] <- gpar(fontsize = 15)
gridfortable.2$grobs[ind11][[1]][["gp"]] <- gpar(fill = "darkolivegreen1", col = "darkolivegreen4", lwd=5)
gridfortable.2$grobs[ind12][[1]][["gp"]] <- gpar(fill = "darkolivegreen1", col = "darkolivegreen4", lwd=5)
gridfortable.2$grobs[ind13][[1]][["gp"]] <- gpar(fill = "darkolivegreen1", col = "darkolivegreen4", lwd=5)
gridfortable.2$grobs[ind14][[1]][["gp"]] <- gpar(fill = "darkolivegreen1", col = "darkolivegreen4", lwd=5)
gridfortable.2$grobs[ind15][[1]][["gp"]] <- gpar(fill = "darkolivegreen1", col = "darkolivegreen4", lwd=5)
gridfortable.2$grobs[ind16][[1]][["gp"]] <- gpar(fill = "darkolivegreen1", col = "darkolivegreen4", lwd=5)
gridfortable.2$grobs[ind17][[1]][["gp"]] <- gpar(fill = "darkolivegreen1", col = "darkolivegreen4", lwd=5)
gridfortable.2$grobs[ind31][[1]][["gp"]] <- gpar(fill = "cadetblue1",      col = "darkolivegreen4", lwd=5)
gridfortable.2$grobs[ind32][[1]][["gp"]] <- gpar(fill = "cadetblue1",      col = "darkolivegreen4", lwd=5)
gridfortable.2$grobs[ind33][[1]][["gp"]] <- gpar(fill = "cadetblue1",      col = "darkolivegreen4", lwd=5)
gridfortable.2$grobs[ind34][[1]][["gp"]] <- gpar(fill = "cadetblue1",      col = "darkolivegreen4", lwd=5)
gridfortable.2$grobs[ind35][[1]][["gp"]] <- gpar(fill = "cadetblue1",      col = "darkolivegreen4", lwd=5)
gridfortable.2$grobs[ind36][[1]][["gp"]] <- gpar(fill = "cadetblue1",      col = "darkolivegreen4", lwd=5)
gridfortable.2$grobs[ind37][[1]][["gp"]] <- gpar(fill = "cadetblue1",      col = "darkolivegreen4", lwd=5)
gridfortable.2$grobs[ind38][[1]][["gp"]] <- gpar(fill = "cadetblue1",      col = "darkolivegreen4", lwd=5)
gridfortable.2$grobs[ind39][[1]][["gp"]] <- gpar(fill = "cadetblue1",      col = "darkolivegreen4", lwd=5)
grid.newpage()
grid.draw(gridfortable.2)

# compute statistics for sum of flights
# frac RadIA processed


