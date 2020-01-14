#----------------------------
#
# Name:    MIXPHA_assess_RadIA.R
# 
# Created: 2.8.2019 dserke
#
#----------------------------

# read in nc file
# start with KGRB volume between 18:48 and 19:03 GMT (3 volumes) on 2/7/2019
radar.tilt         <- nc_open("/d1/serke/projects/case_studies/ICICLE_2018/data/nc_nexrad/20190207_F07/20190207_1902Zvol_KGRB.nc", write = FALSE, verbose = TRUE)
print(paste("The file has", radar.tilt$nvars, "variables"))
radar.tilt.var.num <- c(1:6)
for (i in 1:length(radar.tilt.var.num)) {
  radar.tilt.nam <- paste("v", radar.tilt.var.num[i], sep = "")
  assign(radar.tilt.nam, radar.tilt$var[[radar.tilt.var.num[i]]])
  #print(paste("V", i, " has the name", v1$name))
}

DBZ                      <- ncvar_get(radar.tilt, v1 )
VEL                      <- ncvar_get(radar.tilt, v2 )
SW                       <- ncvar_get(radar.tilt, v3 )
ZDR                      <- ncvar_get(radar.tilt, v4 )
PHI                      <- ncvar_get(radar.tilt, v5 )
RHO                      <- ncvar_get(radar.tilt, v6 )

nc_close(radar.tilt)

# load NRC Convair alt/temp profile
alt.m  <- c( 332,  560,  760, 1034, 1278, 1542, 1652, 1777, 1953, 2273, 2505, 2752, 3010, 3244, 3523, 3754, 4012, 4223, 4493, 4766, 5002)
temp.c <- c(-4.2, -4.2, -5.5, -6.2, -5.6, -2.5,  0.7,  1.1, 1.35,  0.1, -0.7, -1.9, -3.4, -4.6, -6.3, -7.1, -9.0, -10.7, -12.2, -13.6, -15.1)
NRC.df <- data.frame(alt.m, temp.c)

# secondary KMQT 13:52 to 13:54 GMT on 2/4/2019

# convert range index to alt
range.km     <- seq(from = 0.250, to = 458, by = 0.250)
elev.deg     <- 19.7
elev.rad     <- NISTdegTOradian(elev.deg)
alt.km       <- range.km * tan(elev.rad)

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
for (i in 1:dim(DZ)[1]) {
  DBZ.mean.StoW[i] <- mean(DBZ[i, 360:540, 10], na.rm = TRUE)
  DBZ.mean.WtoN[i] <- mean(DBZ[i, 540:720, 10], na.rm = TRUE)
  DBZ.mean.NtoE[i] <- mean(DBZ[i, 1  :180, 10], na.rm = TRUE)
  DBZ.mean.EtoS[i] <- mean(DBZ[i, 180:360, 10], na.rm = TRUE)
  ZDR.mean.StoW[i] <- mean(ZDR[i, 360:540, 10], na.rm = TRUE)
  ZDR.mean.WtoN[i] <- mean(ZDR[i, 540:720, 10], na.rm = TRUE)
  ZDR.mean.NtoE[i] <- mean(ZDR[i, 1  :180, 10], na.rm = TRUE)
  ZDR.mean.EtoS[i] <- mean(ZDR[i, 180:360, 10], na.rm = TRUE)
  RHO.mean.StoW[i] <- mean(RHO[i, 360:540, 10], na.rm = TRUE)
  RHO.mean.WtoN[i] <- mean(RHO[i, 540:720, 10], na.rm = TRUE)
  RHO.mean.NtoE[i] <- mean(RHO[i, 1  :180, 10], na.rm = TRUE)
  RHO.mean.EtoS[i] <- mean(RHO[i, 180:360, 10], na.rm = TRUE)
}

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

# plot for FZRA detect
plot.new()
par(mfrow = c(1, 3))
plot(DBZ.mean.StoW[1:120],  alt.km[1:120], main = "2019/02/07, 19:02 UTC, KGRB, DBZ", type = "l", xlim=c(-20, 40), xlab = "REFL [dBZ]", ylab = "Alt [km]", col = "blue", lwd = 2)
lines(DBZ.mean.WtoN[1:120], alt.km[1:120], col = "red")
lines(DBZ.mean.NtoE[1:120], alt.km[1:120], col = "green")
lines(DBZ.mean.EtoS[1:120], alt.km[1:120], col = "yellow")
lines(NRC.df$temp.c,        NRC.df$alt.m/1000,  col = "black")
#abline(h = 1.35, col = "blue",   lty = 2)
abline(h = 4.00, col = "black",  lty = 2)
abline(h = 5.00, col = "black",  lty = 2)
legend("topright", c("S-W", "W-N", "N-E", "E-S", "AC T [C]"), col=1:4)
grid()
plot(ZDR.mean.StoW[1:120],  alt.km[1:120], main = "2019/02/07, 19:02 UTC, KGRB, ZDR", type = "l", xlab = "ZDR [dB]", ylab = "Alt [km]", col = "blue", lwd = 2)
lines(ZDR.mean.WtoN[1:120], alt.km[1:120], col = "red")
lines(ZDR.mean.NtoE[1:120], alt.km[1:120], col = "green")
lines(ZDR.mean.EtoS[1:120], alt.km[1:120], col = "yellow")
#abline(h = 1.35, col = "blue",   lty = 2)
abline(h = 4.00, col = "black",  lty = 2)
abline(h = 5.00, col = "black",  lty = 2)
legend("topleft", c("S-W", "W-N", "N-E", "E-S"))
grid()
plot(RHO.mean.StoW[1:120],  alt.km[1:120], main = "2019/02/07, 19:02 UTC, KGRB, RHO", type = "l", xlab = "RHO [ ]", ylab = "Alt [km]", col = "blue", lwd = 2)
lines(RHO.mean.WtoN[1:120], alt.km[1:120], col = "red")
lines(RHO.mean.NtoE[1:120], alt.km[1:120], col = "green")
lines(RHO.mean.EtoS[1:120], alt.km[1:120], col = "yellow")
#abline(h = 1.35, col = "blue",   lty = 2)
abline(h = 4.00, col = "black",  lty = 2)
abline(h = 5.00, col = "black",  lty = 2)
legend("topleft", c("S-W", "W-N", "N-E", "E-S"), col=1:4)
grid()

