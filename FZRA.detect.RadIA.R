#----------------------------
#
# Name:    FZRA_detect_RadIA.R
# 
# Created: 2.8.2019 dserke
#
#----------------------------

# read in nc file
# KGRB volume between 18:48 and 19:03 GMT (3 volumes) on 2/7/2019
#radar.vol         <- nc_open("/d1/serke/projects/case_studies/ICICLE_2018/data/nc_nexrad/KGRB/polar/20190207/190217.nc", write = FALSE, verbose = TRUE)
# load NRC Convair alt/temp profile
alt.m      <- c( 332,  560,  760, 1034, 1278, 1542, 1652, 1777, 1953, 2273, 2505, 2752, 3010, 3244, 3523, 3754, 4012, 4223, 4493, 4766, 5002)
temp.c     <- c(-4.2, -4.2, -5.5, -6.2, -5.6, -2.5,  0.7,  1.1, 1.35,  0.1, -0.7, -1.9, -3.4, -4.6, -6.3, -7.1, -9.0, -10.7, -12.2, -13.6, -15.1)
NRC.F07.df <- data.frame(alt.m, temp.c)
# KMQT volume between 13:51 and 13:53 GMT (1 volume) on 2/4/2019
radar.vol         <- nc_open("/d1/serke/projects/case_studies/ICICLE_2018/data/nc_nexrad/KMQT/polar/20190204/134405.nc", write = FALSE, verbose = TRUE)
#radar.vol         <- nc_open("/d1/serke/projects/case_studies/ICICLE_2018/data/nc_nexrad/KMQT/polar/20190204/135012.nc", write = FALSE, verbose = TRUE)
#radar.vol         <- nc_open("/d1/serke/projects/case_studies/ICICLE_2018/data/nc_nexrad/KMQT/polar/20190204/135619.nc", write = FALSE, verbose = TRUE)
# secondary KMQT 13:52 to 13:54 GMT on F04 2/4/2019
alt.m      <- c(6323, 6306, 6159, 6117, 6034, 5978, 5906, 5842, 5781, 5705, 5644, 5496, 5408, 5336, 5236, 5138, 4996, 4884, 4819, 4796, 4743, 4614, 4498, 4408, 4308, 4226, 4155, 4061, 3959, 3903, 3802, 3675, 3627, 3533, 3353, 3280, 3200, 3084, 3055, 2932, 2844, 2791, 2690, 2596, 2526, 2441, 2371, 2304, 2189, 2038, 1990, 1998, 1913, 1838, 1738, 1676, 1660, 1591, 1519, 1466, 1354, 1280, 1233, 1177, 1172, 1176, 1182, 1174, 1180, 1178, 1183, 1175, 1164, 1157, 1134, 1189, 1180, 1191, 1180, 1172, 1148, 1084, 980, 886, 762, 664)
temp.c     <- c(-23.9, -23.9, -22.9, -22.0, -21.5, -21.0, -20.4, -19.9, -20.7, -20.2, -19.9, -19.15, -18.1, -17.6, -16.7, -16.2, -15.3, -14.5, -14.1, -13.8, -13.7, -12.8, -11.7, -11.3, -10.7, -10.1, -9.7, -9.5, -8.9, -8.6, -7.9, -7.6, -7.0, -6.9, -5.6, -5.3, -4.8, -4.1, -4.0, -3.2, -2.4, -2.2, -1.9, -1.3, -1.1, -1.2, -0.4, -0.4, 0.1, 0.8, 1.2, 1.1, 0.9, 1.9, 2.3, 2.3, 2.3, 2.2, 2.5, 2.7, 2.6, 2.1, 1.1, 1.0, 0.8, 1.2, 1.7, 2.1, 1.7, 1.4, 1.5, 1.3, 1.4, 1.7, 0.5, 1.3, 1.0, 1.3, 0.7, 0.5, 0.5, 0.4, -0.1, -2.3, -2.7, -2.5)
#alt.m      <- c(1192, 1180, 1172, 1148, 1084,  980,  886,  762,  664)
#temp.c     <- c( 1.4,  0.8,  0.5,  0.6,  0.4, -0.2, -2.3, -2.7, -2.6)
NRC.F04.df <- data.frame(alt.m, temp.c)

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

# convert range index to alt
range.km     <- seq(from = 0.250, to = 458, by = 0.250)
elev.deg     <- 19.7
elev.rad     <- NISTdegTOradian(elev.deg)
alt.km       <- range.km * tan(elev.rad)

# calc KDP, defined as the moving slope of the least squared fit to PHIDP
# THIS IS NOT WORKING CORRECTLY YET  ######
ls.fit.phi   <- lsfit(x = range.km, y = PHI[, 180, dim(PHI)[3]])
KDP          <- movingFun(ls.fit.phi$residuals, n=5, na.rm=TRUE)
plot.new()
plot(range.km[1:190], KDP[1:190], type="b")

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
KDP.mean.StoW <- rep(0, length(KDP))
KDP.mean.WtoN <- rep(0, length(KDP))
KDP.mean.NtoE <- rep(0, length(KDP))
KDP.mean.EtoS <- rep(0, length(KDP))
for (i in 1:dim(DZ)[1]) {
  DBZ.mean.StoW[i] <- mean(DBZ[i, 360:540, 15], na.rm = TRUE)
  DBZ.mean.WtoN[i] <- mean(DBZ[i, 540:720, 15], na.rm = TRUE)
  DBZ.mean.NtoE[i] <- mean(DBZ[i, 1  :180, 15], na.rm = TRUE)
  DBZ.mean.EtoS[i] <- mean(DBZ[i, 180:360, 15], na.rm = TRUE)
  ZDR.mean.StoW[i] <- mean(ZDR[i, 360:540, 15], na.rm = TRUE)
  ZDR.mean.WtoN[i] <- mean(ZDR[i, 540:720, 15], na.rm = TRUE)
  ZDR.mean.NtoE[i] <- mean(ZDR[i, 1  :180, 15], na.rm = TRUE)
  ZDR.mean.EtoS[i] <- mean(ZDR[i, 180:360, 15], na.rm = TRUE)
  RHO.mean.StoW[i] <- mean(RHO[i, 360:540, 15], na.rm = TRUE)
  RHO.mean.WtoN[i] <- mean(RHO[i, 540:720, 15], na.rm = TRUE)
  RHO.mean.NtoE[i] <- mean(RHO[i, 1  :180, 15], na.rm = TRUE)
  RHO.mean.EtoS[i] <- mean(RHO[i, 180:360, 15], na.rm = TRUE)
  KDP.mean.StoW[i] <- mean(KDP[360:540], na.rm = TRUE)
  KDP.mean.WtoN[i] <- mean(KDP[540:720], na.rm = TRUE)
  KDP.mean.NtoE[i] <- mean(KDP[1  :180], na.rm = TRUE)
  KDP.mean.EtoS[i] <- mean(KDP[180:360], na.rm = TRUE)
}

# enter criteria for detection from VDBetal2016

# estimate FZRA rate
# find 1km alt DBZ, ZDR, RHO, KDP
DBZ.mean.StoW.1500m <- DBZ.mean.StoW[6]
ZDR.mean.StoW.1000m <- ZDR.mean.StoW[4]
RHO.mean.StoW.1000m <- RHO.mean.StoW[4]
KDP.mean.StoW.1000m <- KDP.mean.StoW[4]
# print values
print(paste("mean DBZ S-to-W @ 1.5km alt = ", DBZ.mean.StoW.1500m, sep=""))
print("VDBetal2016 DBZ @ 1.5km alt r = 0.60")
print(paste("mean ZDR S-to-W @ 1.0km alt = ", ZDR.mean.StoW.1000m, sep=""))
print("VDBetal2016 ZDR @ 1.0km alt r = 0.64")
print(paste("mean RHO S-to-W @ 1.0km alt = ", RHO.mean.StoW.1000m, sep=""))
print("VDBetal2016 RHO @ 1.0km alt r = 0.21")
print(paste("mean KDP S-to-W @ 1.0km alt = ", KDP.mean.StoW.1000m, sep=""))
print("VDBetal2016 KDP @ 1.0km alt r = 0.76")
# compare to VDBetal2016 best fit lines, pulled these values off of the plots on Figure 6 of VDBetal2016
m.zhh.1500m <-  2.2
m.kdp.1000m <-  0.0125
m.zdr.1000m <-  0.22
b.zhh.1500m <- 26.0
b.kdp.1000m <- -0.01
b.zdr.1000m <-  0.1
# y=m*x+b, so solve for x ... x=(y-b)/m
FZRA.rate.zhh.1500m <- (DBZ.mean.StoW.1500m - b.zhh.1500m) / m.zhh.1500m
FZRA.rate.zdr.1000m <- (ZDR.mean.StoW.1000m - b.zdr.1000m) / m.zdr.1000m
FZRA.rate.kdp.1000m <- (KDP.mean.StoW.1000m - b.kdp.1000m) / m.kdp.1000m
print(paste("FRZA-rate(DBZ) = ", signif(FZRA.rate.zhh.1500m, 2), " [mm hr-1]", sep=""))
print(paste("FRZA-rate(ZDR) = ", signif(FZRA.rate.zdr.1000m, 2), " [mm hr-1]", sep=""))

# plot for FZRA detect on 20190207 @ 19:02 GMT @ KGRB
plot.new()
par(mfrow = c(1, 3))
plot(DBZ.mean.StoW[1:120],  alt.km[1:120], main = "2019/02/07, 19:02 UTC, KGRB, DBZ", type = "l", xlab = "REFL [dBZ]", ylab = "Alt [km]", col = "blue", lwd = 2)
lines(DBZ.mean.WtoN[1:120], alt.km[1:120], col = "red")
lines(DBZ.mean.NtoE[1:120], alt.km[1:120], col = "green")
lines(DBZ.mean.EtoS[1:120], alt.km[1:120], col = "yellow")
lines(NRC.F07.df$temp.c,        NRC.F07.df$alt.m/1000,  col = "black")
abline(h = 1.35, col = "blue",   lty = 2)
abline(h = 2.28, col = "black",  lty = 2)
abline(h = 1.60, col = "black",  lty = 2)
legend("topright", c("S-W", "W-N", "N-E", "E-S", "AC T [C]"), col=1:4)
grid()
plot(ZDR.mean.StoW[1:120],  alt.km[1:120], main = "2019/02/07, 19:02 UTC, KGRB, ZDR", type = "l", xlab = "ZDR [dB]", ylab = "Alt [km]", col = "blue", lwd = 2)
lines(ZDR.mean.WtoN[1:120], alt.km[1:120], col = "red")
lines(ZDR.mean.NtoE[1:120], alt.km[1:120], col = "green")
lines(ZDR.mean.EtoS[1:120], alt.km[1:120], col = "yellow")
abline(h = 1.35, col = "blue",   lty = 2)
abline(h = 2.28, col = "black",  lty = 2)
abline(h = 1.60, col = "black",  lty = 2)
legend("topleft", c("S-W", "W-N", "N-E", "E-S"))
grid()
plot(RHO.mean.StoW[1:120],  alt.km[1:120], main = "2019/02/07, 19:02 UTC, KGRB, RHO", type = "l", xlab = "RHO [ ]", ylab = "Alt [km]", col = "blue", lwd = 2)
lines(RHO.mean.WtoN[1:120], alt.km[1:120], col = "red")
lines(RHO.mean.NtoE[1:120], alt.km[1:120], col = "green")
lines(RHO.mean.EtoS[1:120], alt.km[1:120], col = "yellow")
abline(h = 1.35, col = "blue",   lty = 2)
abline(h = 2.28, col = "black",  lty = 2)
abline(h = 1.60, col = "black",  lty = 2)
legend("topleft", c("S-W", "W-N", "N-E", "E-S"), col=1:4)
grid()

# plot for FZRA detect on 20190204 @ 13:50 GMT @ KMQT
plot.new()
par(mfrow = c(1, 3))
plot(DBZ.mean.StoW[1:120],  alt.km[1:120], main = "2019/02/04, 13:44 UTC, KMQT, DBZ", type = "l", xlab = "REFL [dBZ]", ylab = "Alt [km]", col = "blue", lwd = 2)
#plot(DBZ.mean.StoW[1:120],  alt.km[1:120], main = "2019/02/04, 13:50 UTC, KMQT, DBZ", xlim=c(-20, 40), type = "l", xlab = "REFL [dBZ]", ylab = "Alt [km]", col = "blue", lwd = 1)
#plot(DBZ.mean.StoW[1:120],  alt.km[1:120], main = "2019/02/04, 13:56 UTC, KMQT, DBZ", type = "l", xlab = "REFL [dBZ]", ylab = "Alt [km]", col = "blue", lwd = 2)
lines(DBZ.mean.WtoN[1:120], alt.km[1:120],          col = "red")
lines(DBZ.mean.NtoE[1:120], alt.km[1:120],          col = "green")
lines(DBZ.mean.EtoS[1:120], alt.km[1:120],          col = "yellow")
lines(NRC.F04.df$temp.c,    NRC.F04.df$alt.m/1000,  col = "black")
abline(h = 3.60, col = "red",    lty = 2)
abline(h = 2.20, col = "black",  lty = 2)
abline(h = 1.00, col = "black",  lty = 2)
abline(h = 0.98, col = "blue",   lty = 2)
legend("topright", c("S-W", "W-N", "N-E", "E-S", "AC T [C]"), col=1:4)
grid()
plot(ZDR.mean.StoW[1:120],  alt.km[1:120], main = "2019/02/04, 13:44 UTC, KMQT, ZDR", type = "l", xlab = "ZDR [dB]", ylab = "Alt [km]", col = "blue", lwd = 2)
#plot(ZDR.mean.StoW[1:120],  alt.km[1:120], main = "2019/02/04, 13:50 UTC, KMQT, ZDR", type = "l", xlab = "ZDR [dB]", ylab = "Alt [km]", col = "blue", lwd = 1)
#plot(ZDR.mean.StoW[1:120],  alt.km[1:120], main = "2019/02/04, 13:56 UTC, KMQT, ZDR", type = "l", xlab = "ZDR [dB]", ylab = "Alt [km]", col = "blue", lwd = 2)
lines(ZDR.mean.WtoN[1:120], alt.km[1:120], col = "red")
lines(ZDR.mean.NtoE[1:120], alt.km[1:120], col = "green")
lines(ZDR.mean.EtoS[1:120], alt.km[1:120], col = "yellow")
abline(h = 4.20, col = "red",    lty = 2)
abline(h = 2.20, col = "black",  lty = 2)
abline(h = 1.00, col = "black",  lty = 2)
abline(h = 0.80, col = "blue",   lty = 2)
legend("topleft", c("S-W", "W-N", "N-E", "E-S"))
grid()
plot(RHO.mean.StoW[1:120],  alt.km[1:120], main = "2019/02/04, 13:44 UTC, KMQT, RHO", type = "l", xlab = "RHO [ ]", ylab = "Alt [km]", col = "blue", lwd = 2)
#plot(RHO.mean.StoW[1:120],  alt.km[1:120], main = "2019/02/04, 13:50 UTC, KMQT, RHO", type = "l", xlab = "RHO [ ]", ylab = "Alt [km]", col = "blue", lwd = 1)
#plot(RHO.mean.StoW[1:120],  alt.km[1:120], main = "2019/02/04, 13:56 UTC, KMQT, RHO", type = "l", xlab = "RHO [ ]", ylab = "Alt [km]", col = "blue", lwd = 2)
lines(RHO.mean.WtoN[1:120], alt.km[1:120], col = "red")
lines(RHO.mean.NtoE[1:120], alt.km[1:120], col = "green")
lines(RHO.mean.EtoS[1:120], alt.km[1:120], col = "yellow")
abline(h = 4.20, col = "red",    lty = 2)
abline(h = 2.20, col = "black",  lty = 2)
abline(h = 1.00, col = "black",  lty = 2)
abline(h = 0.80, col = "blue",   lty = 2)
legend("topleft", c("S-W", "W-N", "N-E", "E-S"), col=1:4)
grid()

