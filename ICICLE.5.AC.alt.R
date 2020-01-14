# info on missed approach/profile during F24

# unit conversion
library(NISTunits)

# define NRC CNV flight track (planet) csv file path
flight.track.planet.dir            <- file.path("/d1/serke/projects/case_studies/ICICLE_2018/data/flight_track_planet/")
flight.track.df                    <- read.csv(paste(flight.track.planet.dir, "F21_NRC_track_2019-02-24.csv", sep = ""), header = FALSE, sep = ",", dec = ".", stringsAsFactors=FALSE)
colnames(flight.track.df)          <- c("type", "ac.date.posix", "lat", "lon", "alt.m", "", "alt.ft", "", "radar.alt.ft", "ground.spd.mps", "true.spd.mps", "ind.spd.kts", "mach.num", "vert.vel.mps", "true.hdg.deg", "track.deg", "drift.deg", "pitch.deg", "roll.deg", "side.slip.deg", "angle.atk.deg", "amb.temp.deg.C", "dewpt.temp.deg.C", "", "", "", "?5", "?6")
head(flight.track.df)

# find and display specific times in the flight track
flight.track.df$ac.date.posix[256:330]

# sample plot
plot(flight.track.df$alt.m[256:330])
grid()

  radia.SRPC.nc.ALL.filename        <- paste(radia.SRPC.nc.ALL.dir, "20190224_123144.nc", sep="")
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
  
  
# calc alt ARL based on theta and range index
range.km    <- seq(from = 1, to = 110, by = 0.25)
alt.ARL.km  <- range.km * sin(NISTdegTOradian(0.5))
alt.ARL.ft  <- alt.ARL.km * 3281
  

plot(FZDZ.SRPC[1:500,1,1], ylim=c(0,1))
