#----------------------------
#
# Name:    ZDRBB_assess_RadIA.R
# 
# Created: 2.8.2019 dserke
#
#----------------------------

library(pspline)

# read in radar moms nc file
# KMQT volume between 13:25 and at 15:02 GMT (2 volumes) on 2/4/2019
radar.vol         <- nc_open("/d1/serke/projects/case_studies/ICICLE_2018/data/nc_nexrad/KMQT/polar/20190204/132543.nc", write = FALSE, verbose = TRUE)
print(paste("The file has", radar.vol$nvars, "variables"))
radar.vol.var.num <- c(1:6)
for (i in 1:length(radar.vol.var.num)) {
  radar.vol.nam <- paste("v", radar.vol.var.num[i], sep = "")
  assign(radar.vol.nam, radar.vol$var[[radar.vol.var.num[i]]])
  #print(paste("V", i, " has the name", v1$name))
}
DBZ                      <- ncvar_get(radar.vol, v1 )
VEL                      <- ncvar_get(radar.vol, v2 )
SW                       <- ncvar_get(radar.vol, v3 )
ZDR                      <- ncvar_get(radar.vol, v4 )
PHI                      <- ncvar_get(radar.vol, v5 )
RHO                      <- ncvar_get(radar.vol, v6 )
nc_close(radar.vol)

# read in radia vol nc file
radia.vol         <- nc_open("/d1/serke/projects/case_studies/ICICLE_2018/data/nc_radia/interest/KMQT/polar/20190204/132543.nc", write = FALSE, verbose = TRUE)
print(paste("The file has", radia.vol$nvars, "variables"))
radia.vol.var.num <- c(1:5)
for (i in 1:length(radia.vol.var.num)) {
  radia.vol.nam <- paste("v", radia.vol.var.num[i], sep = "")
  assign(radia.vol.nam, radia.vol$var[[radia.vol.var.num[i]]])
  #print(paste("V", i, " has the name", v1$name))
}
FRZDRZ                   <- ncvar_get(radia.vol, v1 )
SLW                      <- ncvar_get(radia.vol, v2 )
MIXPHA                   <- ncvar_get(radia.vol, v3 )
PLATES                   <- ncvar_get(radia.vol, v4 )
RADIA2                   <- ncvar_get(radia.vol, v5 )
nc_close(radia.vol)

# load NRC Convair alt/temp profile
# secondary KMQT 13:25 to 13:54 GMT on F04 2/4/2019
alt.m      <- c(6323, 6306, 6159, 6117, 6034, 5978, 5906, 5842, 5781, 5705, 5644, 5496, 5408, 5336, 5236, 5138, 4996, 4884, 4819, 4796, 4743, 4614, 4498, 4408, 4308, 4226, 4155, 4061, 3959, 3903, 3802, 3675, 3627, 3533, 3353, 3280, 3200, 3084, 3055, 2932, 2844, 2791, 2690, 2596, 2526, 2441, 2371, 2304, 2189, 2038, 1990, 1998, 1913, 1838, 1738, 1676, 1660, 1591, 1519, 1466, 1354, 1280, 1233, 1177, 1172, 1176, 1182, 1174, 1180, 1178, 1183, 1175, 1164, 1157, 1134, 1189, 1180, 1191, 1180, 1172, 1148, 1084, 980, 886, 762, 664)
temp.c     <- c(-23.9, -23.9, -22.9, -22.0, -21.5, -21.0, -20.4, -19.9, -20.7, -20.2, -19.9, -19.15, -18.1, -17.6, -16.7, -16.2, -15.3, -14.5, -14.1, -13.8, -13.7, -12.8, -11.7, -11.3, -10.7, -10.1, -9.7, -9.5, -8.9, -8.6, -7.9, -7.6, -7.0, -6.9, -5.6, -5.3, -4.8, -4.1, -4.0, -3.2, -2.4, -2.2, -1.9, -1.3, -1.1, -1.2, -0.4, -0.4, 0.1, 0.8, 1.2, 1.1, 0.9, 1.9, 2.3, 2.3, 2.3, 2.2, 2.5, 2.7, 2.6, 2.1, 1.1, 1.0, 0.8, 1.2, 1.7, 2.1, 1.7, 1.4, 1.5, 1.3, 1.4, 1.7, 0.5, 1.3, 1.0, 1.3, 0.7, 0.5, 0.5, 0.4, -0.1, -2.3, -2.7, -2.5)
#alt.m      <- c(1192, 1180, 1172, 1148, 1084,  980,  886,  762,  664)
#temp.c     <- c( 1.4,  0.8,  0.5,  0.6,  0.4, -0.2, -2.3, -2.7, -2.6)
NRC.F04.df <- data.frame(alt.m, temp.c)

# convert nexrad moment range index to alt
range.km       <- seq(from = 0.250, to = 458, by = 0.250)
elev.radar.deg <- 19.7
elev.radar.rad <- NISTdegTOradian(elev.radar.deg)
alt.radar.km   <- range.km * tan(elev.radar.rad)

# convert radia moment range index to alt
range.km       <- seq(from = 0.250, to = 458, by = 0.250)
elev.radia.deg <- 6.4
elev.radia.rad <- NISTdegTOradian(elev.radia.deg)
alt.radia.km   <- range.km * tan(elev.radia.rad)

# calc KDP here
# KDP is slope of least squared fit of PHIDP
#ls.fit.phi <- lsfit(x=range.km, y=PHI[,180])

# quality control the data
DBZ[DBZ < -50]             <- NA
ZDR[ZDR < -5 & ZDR > 9]    <- NA
RHO[RHO < 0.5 & RHO > 1.1] <- NA

# mean data at each range, which is a QVP
DBZ.mean.StoW <- rep(0, dim(DBZ)[1])
DBZ.mean.WtoN <- rep(0, dim(DBZ)[1])
DBZ.mean.NtoE <- rep(0, dim(DBZ)[1])
DBZ.mean.EtoS <- rep(0, dim(DBZ)[1])
ZDR.mean.StoW <- rep(0, dim(ZDR)[1])
ZDR.mean.WtoN <- rep(0, dim(ZDR)[1])
ZDR.mean.NtoE <- rep(0, dim(ZDR)[1])
ZDR.mean.EtoS <- rep(0, dim(ZDR)[1])
RHO.mean.StoW <- rep(0, dim(RHO)[1])
RHO.mean.WtoN <- rep(0, dim(RHO)[1])
RHO.mean.NtoE <- rep(0, dim(RHO)[1])
RHO.mean.EtoS <- rep(0, dim(RHO)[1])
PHI.mean.StoW <- rep(0, dim(PHI)[1])
PHI.mean.WtoN <- rep(0, dim(PHI)[1])
PHI.mean.NtoE <- rep(0, dim(PHI)[1])
PHI.mean.EtoS <- rep(0, dim(PHI)[1])
MIXPHA.mean.StoW <- rep(0, dim(MIXPHA)[1])
MIXPHA.mean.WtoN <- rep(0, dim(MIXPHA)[1])
MIXPHA.mean.NtoE <- rep(0, dim(MIXPHA)[1])
MIXPHA.mean.EtoS <- rep(0, dim(MIXPHA)[1])
FRZDRZ.mean.StoW <- rep(0, dim(FRZDRZ)[1])
FRZDRZ.mean.WtoN <- rep(0, dim(FRZDRZ)[1])
FRZDRZ.mean.NtoE <- rep(0, dim(FRZDRZ)[1])
FRZDRZ.mean.EtoS <- rep(0, dim(FRZDRZ)[1])
PLATES.mean.StoW <- rep(0, dim(PLATES)[1])
PLATES.mean.WtoN <- rep(0, dim(PLATES)[1])
PLATES.mean.NtoE <- rep(0, dim(PLATES)[1])
PLATES.mean.EtoS <- rep(0, dim(PLATES)[1])
SLW.mean.StoW <- rep(0, dim(SLW)[1])
SLW.mean.WtoN <- rep(0, dim(SLW)[1])
SLW.mean.NtoE <- rep(0, dim(SLW)[1])
SLW.mean.EtoS <- rep(0, dim(SLW)[1])

tilt.num.radar   <- dim(DBZ)[3]
tilt.num.radia   <- dim(FRZDRZ)[3]

for (i in 1:dim(DBZ)[1]) {
  DBZ.mean.StoW[i]    <- mean(DBZ[i, 360:540, tilt.num.radar], na.rm = TRUE)
  DBZ.mean.WtoN[i]    <- mean(DBZ[i, 540:720, tilt.num.radar], na.rm = TRUE)
  DBZ.mean.NtoE[i]    <- mean(DBZ[i, 1  :180, tilt.num.radar], na.rm = TRUE)
  DBZ.mean.EtoS[i]    <- mean(DBZ[i, 180:360, tilt.num.radar], na.rm = TRUE)
  ZDR.mean.StoW[i]    <- mean(ZDR[i, 360:540, tilt.num.radar], na.rm = TRUE)
  ZDR.mean.WtoN[i]    <- mean(ZDR[i, 540:720, tilt.num.radar], na.rm = TRUE)
  ZDR.mean.NtoE[i]    <- mean(ZDR[i, 1  :180, tilt.num.radar], na.rm = TRUE)
  ZDR.mean.EtoS[i]    <- mean(ZDR[i, 180:360, tilt.num.radar], na.rm = TRUE)
  RHO.mean.StoW[i]    <- mean(RHO[i, 360:540, tilt.num.radar], na.rm = TRUE)
  RHO.mean.WtoN[i]    <- mean(RHO[i, 540:720, tilt.num.radar], na.rm = TRUE)
  RHO.mean.NtoE[i]    <- mean(RHO[i, 1  :180, tilt.num.radar], na.rm = TRUE)
  RHO.mean.EtoS[i]    <- mean(RHO[i, 180:360, tilt.num.radar], na.rm = TRUE)
  PHI.mean.StoW[i]    <- mean(PHI[i, 360:540, tilt.num.radar], na.rm = TRUE)
  PHI.mean.WtoN[i]    <- mean(PHI[i, 540:720, tilt.num.radar], na.rm = TRUE)
  PHI.mean.NtoE[i]    <- mean(PHI[i, 1  :180, tilt.num.radar], na.rm = TRUE)
  PHI.mean.EtoS[i]    <- mean(PHI[i, 180:360, tilt.num.radar], na.rm = TRUE)
  MIXPHA.mean.StoW[i] <- mean(MIXPHA[i, 360:540, tilt.num.radia], na.rm = TRUE)
  MIXPHA.mean.WtoN[i] <- mean(MIXPHA[i, 540:720, tilt.num.radia], na.rm = TRUE)
  MIXPHA.mean.NtoE[i] <- mean(MIXPHA[i, 1  :180, tilt.num.radia], na.rm = TRUE)
  MIXPHA.mean.EtoS[i] <- mean(MIXPHA[i, 180:360, tilt.num.radia], na.rm = TRUE)
  FRZDRZ.mean.StoW[i] <- mean(FRZDRZ[i, 360:540, tilt.num.radia], na.rm = TRUE)
  FRZDRZ.mean.WtoN[i] <- mean(FRZDRZ[i, 540:720, tilt.num.radia], na.rm = TRUE)
  FRZDRZ.mean.NtoE[i] <- mean(FRZDRZ[i, 1  :180, tilt.num.radia], na.rm = TRUE)
  FRZDRZ.mean.EtoS[i] <- mean(FRZDRZ[i, 180:360, tilt.num.radia], na.rm = TRUE)
  PLATES.mean.StoW[i] <- mean(PLATES[i, 360:540, tilt.num.radia], na.rm = TRUE)
  PLATES.mean.WtoN[i] <- mean(PLATES[i, 540:720, tilt.num.radia], na.rm = TRUE)
  PLATES.mean.NtoE[i] <- mean(PLATES[i, 1  :180, tilt.num.radia], na.rm = TRUE)
  PLATES.mean.EtoS[i] <- mean(PLATES[i, 180:360, tilt.num.radia], na.rm = TRUE)
  SLW.mean.StoW[i]    <- mean(SLW[i, 360:540, tilt.num.radia], na.rm = TRUE)
  SLW.mean.WtoN[i]    <- mean(SLW[i, 540:720, tilt.num.radia], na.rm = TRUE)
  SLW.mean.NtoE[i]    <- mean(SLW[i, 1  :180, tilt.num.radia], na.rm = TRUE)
  SLW.mean.EtoS[i]    <- mean(SLW[i, 180:360, tilt.num.radia], na.rm = TRUE)
}

# smooth data
DBZ.mean.StoW[is.na(DBZ.mean.StoW)]  <- 0
DBZ.mean.StoW[is.nan(DBZ.mean.StoW)] <- 0
smoothed.DBZ.mean.StoW               <- smooth.Pspline(x = 1:length(DBZ.mean.StoW), y= DBZ.mean.StoW, df=5, method=3)
smoothed.DBZ.mean.StoW.f0            <- predict(smoothed.DBZ.mean.StoW, x=1:length(DBZ.mean.StoW), nderiv=0)
smoothed.DBZ.mean.StoW.f1            <- predict(smoothed.DBZ.mean.StoW, x=1:length(DBZ.mean.StoW), nderiv=1)
smoothed.DBZ.mean.StoW.f2            <- predict(smoothed.DBZ.mean.StoW, x=1:length(DBZ.mean.StoW), nderiv=2)
ZDR.mean.StoW[is.na(ZDR.mean.StoW)]  <- 0
ZDR.mean.StoW[is.nan(ZDR.mean.StoW)] <- 0
smoothed.ZDR.mean.StoW               <- smooth.Pspline(x = 1:length(ZDR.mean.StoW), y= ZDR.mean.StoW, df=5, method=3)
smoothed.ZDR.mean.StoW.f0            <- predict(smoothed.ZDR.mean.StoW, x=1:length(ZDR.mean.StoW), nderiv=0)
smoothed.ZDR.mean.StoW.f1            <- predict(smoothed.ZDR.mean.StoW, x=1:length(ZDR.mean.StoW), nderiv=1)
smoothed.ZDR.mean.StoW.f2            <- predict(smoothed.ZDR.mean.StoW, x=1:length(ZDR.mean.StoW), nderiv=2)

# dendritic growth zone is where first deriv of DBZ is min (slope large neg) and second deriv of ZDR is lt zero
ind.dend.growth.zone                 <- which(smoothed.DBZ.mean.StoW.f1<= -1.5 & smoothed.ZDR.mean.StoW.f2 < 0)

# enter criteria for detection from VDBetal2016

# estimate FZRA rate
# find 1km alt DBZ, ZDR, RHO
DBZ.mean.StoW.1500m <- DBZ.mean.StoW[6]
ZDR.mean.StoW.1000m <- ZDR.mean.StoW[4]
RHO.mean.StoW.1000m <- RHO.mean.StoW[4]
# print values
print(paste("mean DBZ S-to-W @ 1.5km alt = ", DBZ.mean.StoW.1500m, sep=""))
print("VDBetal2016 DBZ @ 1.5km alt r = 0.60")
print(paste("mean ZDR S-to-W @ 1.0km alt = ", ZDR.mean.StoW.1000m, sep=""))
print("VDBetal2016 ZDR @ 1.0km alt r = 0.64")
print(paste("mean RHO S-to-W @ 1.0km alt = ", RHO.mean.StoW.1000m, sep=""))
print("VDBetal2016 RHO @ 1.0km alt r = 0.21")
print("VDBetal2016 KDP @ 1.0km alt r = 0.76")
# compare to VDBetal2016 best fit lines, pulled these values off of the plots on Figure 6 of VDBetal2016
m.zhh.1500m <-  2.2
m.kdp.1000m <-  0.0125
m.zdr.1000m <-  0.22
b.zhh.1500m <- 26.0
b.kdp.1000m <- -0.01
b.zdr.1000m <-  0.1
# y=m*x+b, so x=(y-b)/m
FZRA.rate.zhh.1500m <- (DBZ.mean.StoW.1500m - b.zhh.1500m) / m.zhh.1500m
FZRA.rate.zdr.1000m <- (ZDR.mean.StoW.1000m - b.zdr.1000m) / m.zdr.1000m
print(paste("FRZA-rate(DBZ) = ", signif(FZRA.rate.zhh.1500m, 2), " [mm hr-1]", sep=""))
print(paste("FRZA-rate(ZDR) = ", signif(FZRA.rate.zdr.1000m, 2), " [mm hr-1]", sep=""))

# plot smoothed.dbz, first deriv, second deriv
plot.new()
par(mfrow = c(1, 4))
plot(DBZ.mean.StoW[1:120],             alt.radar.km[1:120], main = "2019/02/04, 13:25 UTC, KMQT, DBZ", type = "l", xlim=c(-20, 40), xlab = "REFL [dBZ]", ylab = "Alt [km]", col = "blue", lwd = 2)
lines(smoothed.DBZ.mean.StoW$y[1:120], alt.radar.km[1:120], col = "blue", lwd=3)
grid()
plot(smoothed.DBZ.mean.StoW.f1[1:120], alt.radar.km[1:120], type="l")
grid()
plot(ZDR.mean.StoW[1:120],             alt.radar.km[1:120], main = "2019/02/04, 13:25 UTC, KMQT, ZDR", type = "l", xlim=c(-1, 6), xlab = "ZDR [dB]", ylab = "Alt [km]", col = "blue", lwd = 2)
grid()
plot(smoothed.ZDR.mean.StoW.f2[1:120], alt.radar.km[1:120], type="l")
grid()
#plot(f2[1:120],alt.radar.km[1:120], type="l")
#grid()

# plot for ZDR BB detect detect
plot.new()
par(mfrow = c(1, 4))
plot(DBZ.mean.StoW[1:120],  alt.radar.km[1:120], main = "2019/02/04, 13:25 UTC, KMQT, DBZ", type = "l", xlim=c(-20, 40), xlab = "REFL [dBZ]", ylab = "Alt [km]", col = "blue", lwd = 2)
lines(smoothed.DBZ.mean.StoW$y[1:120], alt.radar.km[1:120], col = "blue", lwd=3)
lines(DBZ.mean.WtoN[1:120], alt.radar.km[1:120], col = "red")
lines(DBZ.mean.NtoE[1:120], alt.radar.km[1:120], col = "green")
lines(DBZ.mean.EtoS[1:120], alt.radar.km[1:120], col = "yellow")
lines(NRC.df$temp.c,        NRC.df$alt.m/1000,  col = "black")
abline(h = 4.20, col = "blue",   lty = 2)
abline(h = 2.40, col = "black",  lty = 2)
abline(h = 1.70, col = "black",  lty = 2)
legend("topright", c("S-W", "W-N", "N-E", "E-S", "AC T [C]"), col=1:4)
grid()
plot(ZDR.mean.StoW[1:120],  alt.radar.km[1:120], main = "2019/02/04, 13:25 UTC, KMQT, ZDR", type = "l", xlim=c(-1,6), xlab = "ZDR [dB]", ylab = "Alt [km]", col = "blue", lwd = 2)
lines(ZDR.mean.WtoN[1:120], alt.radar.km[1:120], col = "red")
lines(ZDR.mean.NtoE[1:120], alt.radar.km[1:120], col = "green")
lines(ZDR.mean.EtoS[1:120], alt.radar.km[1:120], col = "yellow")
abline(h = 4.20, col = "blue",   lty = 2)
#abline(h = 4.00, col = "black",  lty = 2)
#abline(h = 5.00, col = "black",  lty = 2)
legend("topleft", c("S-W", "W-N", "N-E", "E-S"))
grid()
plot(RHO.mean.StoW[1:120],  alt.radar.km[1:120], main = "2019/02/04, 13:25 UTC, KMQT, RHO", type = "l", xlab = "RHO [ ]", ylab = "Alt [km]", col = "blue", lwd = 2)
lines(RHO.mean.WtoN[1:120], alt.radar.km[1:120], col = "red")
lines(RHO.mean.NtoE[1:120], alt.radar.km[1:120], col = "green")
lines(RHO.mean.EtoS[1:120], alt.radar.km[1:120], col = "yellow")
abline(h = 4.20, col = "blue",   lty = 2)
#abline(h = 4.00, col = "black",  lty = 2)
#abline(h = 5.00, col = "black",  lty = 2)
legend("topleft", c("S-W", "W-N", "N-E", "E-S"), col=1:4)
grid()
plot(PHI.mean.StoW[1:120],  alt.radar.km[1:120], main = "2019/02/04, 13:25 UTC, KMQT, PHI", type = "l", xlab = "PHI [deg]", ylab = "Alt [km]", col = "blue", lwd = 2)
lines(PHI.mean.WtoN[1:120], alt.radar.km[1:120], col = "red")
lines(PHI.mean.NtoE[1:120], alt.radar.km[1:120], col = "green")
lines(PHI.mean.EtoS[1:120], alt.radar.km[1:120], col = "yellow")
abline(h = 4.20, col = "blue",   lty = 2)
#abline(h = 4.00, col = "black",  lty = 2)
#abline(h = 5.00, col = "black",  lty = 2)
legend("topleft", c("S-W", "W-N", "N-E", "E-S"), col=1:4)
grid()

# plot for RadIA Interests
plot.new()
par(mfrow = c(1, 4))
plot(FRZDRZ.mean.StoW[1:250],  alt.radia.km[1:250], main = "2019/02/04, 13:25 UTC, KMQT, FRZDRZ", type = "l", ylim =c(0,11), xlim=c(0, 1), xlab = "INT [0-1]", ylab = "Alt [km]", col = "blue", lwd = 2)
lines(FRZDRZ.mean.WtoN[1:250], alt.radia.km[1:250], col = "red")
lines(FRZDRZ.mean.NtoE[1:250], alt.radia.km[1:250], col = "green")
lines(FRZDRZ.mean.EtoS[1:250], alt.radia.km[1:250], col = "yellow")
#lines(NRC.df$temp.c,        NRC.df$alt.m/1000,  col = "black")
abline(h = 4.20, col = "blue",   lty = 2)
abline(h = 2.40, col = "black",  lty = 2)
abline(h = 1.70, col = "black",  lty = 2)
abline(v = 0.60, col = "black",  lty = 1)
legend("topright", c("S-W", "W-N", "N-E", "E-S"), col=1:4)
grid()
plot(MIXPHA.mean.StoW[1:250],  alt.radia.km[1:250], main = "MIXPHA", type = "l", ylim =c(0,11), xlim=c(0, 1), xlab = "INT [0-1]", ylab = "Alt [km]", col = "blue", lwd = 2)
lines(MIXPHA.mean.WtoN[1:250], alt.radia.km[1:250], col = "red")
lines(MIXPHA.mean.NtoE[1:250], alt.radia.km[1:250], col = "green")
lines(MIXPHA.mean.EtoS[1:250], alt.radia.km[1:250], col = "yellow")
#lines(NRC.df$temp.c,        NRC.df$alt.m/1000,  col = "black")
abline(h = 4.20, col = "blue",   lty = 2)
abline(h = 2.40, col = "black",  lty = 2)
abline(h = 1.70, col = "black",  lty = 2)
abline(v = 0.60, col = "black",  lty = 1)
legend("topright", c("S-W", "W-N", "N-E", "E-S"), col=1:4)
grid()
plot(SLW.mean.StoW[1:250],  alt.radia.km[1:250], main = "smalldrop", type = "l", ylim =c(0,11), xlim=c(0, 1), xlab = "INT [0-1]", ylab = "Alt [km]", col = "blue", lwd = 2)
lines(SLW.mean.WtoN[1:250], alt.radia.km[1:250], col = "red")
lines(SLW.mean.NtoE[1:250], alt.radia.km[1:250], col = "green")
lines(SLW.mean.EtoS[1:250], alt.radia.km[1:250], col = "yellow")
#lines(NRC.df$temp.c,        NRC.df$alt.m/1000,  col = "black")
abline(h = 4.20, col = "blue",   lty = 2)
abline(h = 2.40, col = "black",  lty = 2)
abline(h = 1.70, col = "black",  lty = 2)
abline(v = 0.60, col = "black",  lty = 1)
legend("topright", c("S-W", "W-N", "N-E", "E-S"), col=1:4)
grid()
plot(PLATES.mean.StoW[1:250],  alt.radia.km[1:250], main = "anisotropic", type = "l", ylim =c(0,11), xlim=c(0, 1), xlab = "INT [0-1]", ylab = "Alt [km]", col = "blue", lwd = 2)
lines(PLATES.mean.WtoN[1:250], alt.radia.km[1:250], col = "red")
lines(PLATES.mean.NtoE[1:250], alt.radia.km[1:250], col = "green")
lines(PLATES.mean.EtoS[1:250], alt.radia.km[1:250], col = "yellow")
#lines(NRC.df$temp.c,        NRC.df$alt.m/1000,  col = "black")
abline(h = 4.20, col = "blue",   lty = 2)
abline(h = 2.40, col = "black",  lty = 2)
abline(h = 1.70, col = "black",  lty = 2)
abline(v = 0.80, col = "black",  lty = 1)
legend("topright", c("S-W", "W-N", "N-E", "E-S"), col=1:4)
grid()
