#----------------------------
#
# Name:    SEVERE.SLD.DETECT.R
# 
# Created: 2.24.2019 dserke
#
#----------------------------

library(NISTunits)

# read in nc file
# KIWX volume between 1328 and 1335 GMT (2 volumes) on 2/24/2019
radar.vol.VCP212 <- nc_open("/d1/serke/projects/case_studies/ICICLE_2018/data/nc_nexrad/KIWX/polar/20190224/132817.nc", write = FALSE, verbose = TRUE)
radar.vol.VCP35  <- nc_open("/d1/serke/projects/case_studies/ICICLE_2018/data/nc_nexrad/KIWX/polar/20190224/133552.nc", write = FALSE, verbose = TRUE)

# load NRC Convair alt/temp profile
#alt.m      <- c( 332,  560,  760, 1034, 1278, 1542, 1652, 1777, 1953, 2273, 2505, 2752, 3010, 3244, 3523, 3754, 4012, 4223, 4493, 4766, 5002)
#temp.c     <- c(-4.2, -4.2, -5.5, -6.2, -5.6, -2.5,  0.7,  1.1, 1.35,  0.1, -0.7, -1.9, -3.4, -4.6, -6.3, -7.1, -9.0, -10.7, -12.2, -13.6, -15.1)
#NRC.F07.df <- data.frame(alt.m, temp.c)

print(paste("The file has", radar.vol.VCP35$nvars, "variables"))
radar.vol.var.num <- c(1:6)
for (i in 1:length(radar.vol.var.num)) {
  radar.vol.nam <- paste("v", radar.vol.var.num[i], sep = "")
  assign(radar.vol.nam, radar.vol.VCP35$var[[radar.vol.var.num[i]]])
  #print(paste("V", i, " has the name", v1$name))
}
DBZ.VCP35                      <- ncvar_get(radar.vol.VCP35, v1 )
VEL.VCP35                      <- ncvar_get(radar.vol.VCP35, v2 )
SW.VCP35                       <- ncvar_get(radar.vol.VCP35, v3 )
ZDR.VCP35                      <- ncvar_get(radar.vol.VCP35, v4 )
PHI.VCP35                      <- ncvar_get(radar.vol.VCP35, v5 )
RHO.VCP35                      <- ncvar_get(radar.vol.VCP35, v6 )
nc_close(radar.vol.VCP35)

print(paste("The file has", radar.vol.VCP212$nvars, "variables"))
radar.vol.var.num <- c(1:6)
for (i in 1:length(radar.vol.var.num)) {
  radar.vol.nam <- paste("v", radar.vol.var.num[i], sep = "")
  assign(radar.vol.nam, radar.vol.VCP212$var[[radar.vol.var.num[i]]])
  #print(paste("V", i, " has the name", v1$name))
}
DBZ.VCP212                      <- ncvar_get(radar.vol.VCP212, v1 )
VEL.VCP212                      <- ncvar_get(radar.vol.VCP212, v2 )
SW.VCP212                       <- ncvar_get(radar.vol.VCP212, v3 )
ZDR.VCP212                      <- ncvar_get(radar.vol.VCP212, v4 )
PHI.VCP212                      <- ncvar_get(radar.vol.VCP212, v5 )
RHO.VCP212                      <- ncvar_get(radar.vol.VCP212, v6 )
nc_close(radar.vol.VCP212)

# convert range index to alt
range.km     <- seq(from = 0.575, to = 458.5, by = 0.250)
elev.deg     <- 19.51
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

# create properly sized empty arrays
DBZ.VCP35.std       <- rep(0, dim(DBZ.VCP35)[1])
DBZ.VCP35.mean      <- rep(0, dim(DBZ.VCP35)[1])
DBZ.VCP212.std      <- rep(0, dim(DBZ.VCP212)[1])
DBZ.VCP212.mean     <- rep(0, dim(DBZ.VCP212)[1])
ZDR.VCP35.std       <- rep(0, dim(ZDR.VCP35)[1])
ZDR.VCP35.mean      <- rep(0, dim(ZDR.VCP35)[1])
ZDR.VCP212.std      <- rep(0, dim(ZDR.VCP212)[1])
ZDR.VCP212.mean     <- rep(0, dim(ZDR.VCP212)[1])
RHO.VCP35.std       <- rep(0, dim(RHO.VCP35)[1])
RHO.VCP35.mean      <- rep(0, dim(RHO.VCP35)[1])
RHO.VCP212.std      <- rep(0, dim(RHO.VCP212)[1])
RHO.VCP212.mean     <- rep(0, dim(RHO.VCP212)[1])

# mean data at each range for quadrants, which is a QVP
for (i in 1:dim(DBZ.VCP35)[1]) {
  DBZ.VCP35.std[i]       <- sd(  DBZ.VCP35[i,        , 10], na.rm = TRUE)
  DBZ.VCP35.mean[i]      <- mean(DBZ.VCP35[i,        , 10], na.rm = TRUE)
  ZDR.VCP35.std[i]       <- sd(  ZDR.VCP35[i,        , 10], na.rm = TRUE)
  ZDR.VCP35.mean[i]      <- mean(ZDR.VCP35[i,        , 10], na.rm = TRUE)
  RHO.VCP35.std[i]       <- sd(  RHO.VCP35[i,        , 10], na.rm = TRUE)
  RHO.VCP35.mean[i]      <- mean(RHO.VCP35[i,        , 10], na.rm = TRUE)
}
for (i in 1:dim(DBZ.VCP212)[1]) {
  DBZ.VCP212.std[i]       <- sd(  DBZ.VCP212[i,        , 10], na.rm = TRUE)
  DBZ.VCP212.mean[i]      <- mean(DBZ.VCP212[i,        , 10], na.rm = TRUE)
  ZDR.VCP212.std[i]       <- sd(  ZDR.VCP212[i,        , 10], na.rm = TRUE)
  ZDR.VCP212.mean[i]      <- mean(ZDR.VCP212[i,        , 10], na.rm = TRUE)
  RHO.VCP212.std[i]       <- sd(  RHO.VCP212[i,        , 10], na.rm = TRUE)
  RHO.VCP212.mean[i]      <- mean(RHO.VCP212[i,        , 10], na.rm = TRUE)
}

# plot 
plot.new()
par(mfrow = c(1, 3))
plot(DBZ.VCP35.mean[1:57],      alt.km[1:57], main = "2019/02/24, 13:35 UTC, KIWX, DBZ", type = "l", xlab = "REFL [dBZ]", ylab = "Alt [km]", col = "blue", lwd = 2)
#plot(DBZ.mean.StoW[1:57],  alt.km[1:57], main = "2019/02/24, 13:35 UTC, KIWX, DBZ", type = "l", xlab = "REFL [dBZ]", ylab = "Alt [km]", col = "blue", lwd = 2)
lines(DBZ.VCP212.mean[1:57], alt.km[1:57], col = "green")
#lines(DBZ.mean.NtoE[1:57], alt.km[1:57], col = "green")
#lines(DBZ.mean.EtoS[1:57], alt.km[1:57], col = "yellow")
#lines(NRC.F07.df$temp.c,        NRC.F07.df$alt.m/1000,  col = "black")
#abline(h = 1.85, col = "black",   lty = 2)
abline(h = 1.30, col = "red",     lty = 2)
abline(h = 2.50, col = "red",     lty = 2)
abline(h = 0.00, col = "blue",    lty = 2)
abline(h = 0.92, col = "blue",    lty = 2)
#text(-18, 1.77, "MAX SIG ALT",    col = "black")
text(-18, 1.75, "MAX LWC ALT",    col = "red")
text(-18, 1.38, "MIN LIQ ALT",    col = "red")
text(-18, 2.42, "MAX LIQ ALT",    col = "red")
text(-17, 0.08, "MIN ALLICE ALT", col = "blue")
text(-17, 0.82, "MAX ALLICE ALT", col = "blue")
#legend("topright", c("S-W", "W-N", "N-E", "E-S", "AC T [C]"), col=1:4)
grid()

plot(ZDR.VCP35.mean[1:57],      alt.km[1:57], main = "2019/02/24, 13:35 UTC, KIWX, ZDR", type = "l", xlab = "ZDR [dB]", ylab = "Alt [km]", col = "blue", lwd = 2)
#plot(ZDR.mean.StoW[1:57],  alt.km[1:57], main = "2019/02/24, 13:35 UTC, KIWX, ZDR", type = "l", xlab = "ZDR [dB]", ylab = "Alt [km]", col = "blue", lwd = 2)
lines(ZDR.VCP212.mean[1:57], alt.km[1:57], col = "green")
#lines(ZDR.mean.NtoE[1:57], alt.km[1:57], col = "green")
#lines(ZDR.mean.EtoS[1:57], alt.km[1:57], col = "yellow")
#abline(h = 1.85, col = "black",   lty = 2)
#legend("topleft", c("S-W", "W-N", "N-E", "E-S"))
grid()

plot(RHO.VCP35.mean[1:57],   alt.km[1:57], main = "2019/02/24, 13:35 UTC, KIWX, RHO", type = "l", xlab = "RHO [ ]", ylab = "Alt [km]", col = "blue", lwd = 2)
#plot(RHO.mean.StoW[1:57],  alt.km[1:57], main = "2019/02/24, 13:35 UTC, KIWX, RHO", type = "l", xlab = "RHO [ ]", ylab = "Alt [km]", col = "blue", lwd = 2)
lines(RHO.VCP212.mean[1:57], alt.km[1:57], col = "green")
#lines(RHO.mean.NtoE[1:57], alt.km[1:57], col = "green")
#lines(RHO.mean.EtoS[1:57], alt.km[1:57], col = "yellow")
#abline(h = 1.85, col = "black",   lty = 2)
abline(h = 0.92, col = "blue",    lty = 2)
abline(h = 1.38, col = "red",     lty = 2)
#legend("topleft", c("S-W", "W-N", "N-E", "E-S"), col=1:4)
grid()
