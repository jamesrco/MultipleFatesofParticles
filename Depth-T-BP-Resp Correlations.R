# load in data

Depth.T.BP.Resp.Correlations <- read.csv('/Users/jrcollins/Dropbox/KN207-3 NA VICE/Data/2012 North Atlantic resp & flux manny/Depth-T-BP-Resp Correlations.csv')

# load Type II regression package

library(lmodel2)
library(Hmisc)

# specify symbols and colors for respiration
Sta_Resp <- factor(Depth.T.BP.Resp.Correlations$Sta_Resp)
pch.Resp <- c(25, 25, 25, 25, 25, 25, 23, 23, 23, 23, 23, 23, 24, 24, 24, 24, 24, 24, 24, 24, 21, 22, 14, 24, 23)
col.Resp <- c(1)
bg.Resp <- c(rep(rgb(237,34,36,maxColorValue=255),6),
             rep(rgb(249,191,191,maxColorValue=255),6),
             rep(rgb(153,27,30,maxColorValue=255),8),
             rgb(80,108,179,maxColorValue=255),
             rgb(194,206,233,maxColorValue=255),
             rgb(100,14,14,maxColorValue=255),             
             rgb(153,27,30,maxColorValue=255),
             rgb(249,191,191,maxColorValue=255)
             )

# resp on temperature

# ***** excluding from the regression one observation (16; 1.627410,.5853) which is an obvious outlier

R.on.T.slr <- lm(log(Depth.T.BP.Resp.Correlations$Resp_mg_C_m3_d[-16]) ~ log(Depth.T.BP.Resp.Correlations$Temp_Resp[-16]), weights = log(Depth.T.BP.Resp.Correlations$Uncert_resp_mg_C_m3_d[-16]))
summary(R.on.T.slr)

errbar(y = log(Depth.T.BP.Resp.Correlations$Resp_mg_C_m3_d), x = log(Depth.T.BP.Resp.Correlations$Temp_Resp), col=col.Resp, pch=pch.Resp, bg=bg.Resp, ylab="log (Water col. respiration rate)", xlab="log (Temperature)", yplus = log(Depth.T.BP.Resp.Correlations$Resp_mg_C_m3_d)+log(Depth.T.BP.Resp.Correlations$Uncert_resp_mg_C_m3_d), yminus = log(Depth.T.BP.Resp.Correlations$Resp_mg_C_m3_d)-log(Depth.T.BP.Resp.Correlations$Uncert_resp_mg_C_m3_d))
abline(R.on.T.slr)

pdf(file = "R.on.T_noerrbar.pdf",
    width = 2, height = 2, pointsize = 9,
    bg = "white")
par(pty = "s", mar=c(4,4,1,1))
plot(log(Depth.T.BP.Resp.Correlations$Resp_mg_C_m3_d) ~ log(Depth.T.BP.Resp.Correlations$Temp_Resp), col=col.Resp, pch=pch.Resp, bg=bg.Resp, ylab="log (Water col. respiration rate)", xlab="log (Temperature)",xaxt="n",yaxt="n")
axis(1, tck=-0.05)
axis(2, tck=-0.05) 
abline(R.on.T.slr)
dev.off()

pdf(file = "R.on.T.pdf",
    width = 2, height = 2, pointsize = 9,
    bg = "white")
par(pty = "s", mar=c(4,4,1,1))
errbar(y = log(Depth.T.BP.Resp.Correlations$Resp_mg_C_m3_d), x = log(Depth.T.BP.Resp.Correlations$Temp_Resp), col=col.Resp, pch=pch.Resp, bg=bg.Resp, ylab="log (Water col. respiration rate)", xlab="log (Temperature)", yplus = log(Depth.T.BP.Resp.Correlations$Resp_mg_C_m3_d)+log(Depth.T.BP.Resp.Correlations$Uncert_resp_mg_C_m3_d), yminus = log(Depth.T.BP.Resp.Correlations$Resp_mg_C_m3_d)-log(Depth.T.BP.Resp.Correlations$Uncert_resp_mg_C_m3_d),xaxt="n",yaxt="n")
axis(1, tck=-0.05)
axis(2, tck=-0.05) 
abline(R.on.T.slr)
dev.off()


# resp on depth

R.on.D.slr <- lm(Depth.T.BP.Resp.Correlations$Resp_mg_C_m3_d ~ log(Depth.T.BP.Resp.Correlations$Depth_Resp), weights = Depth.T.BP.Resp.Correlations$Uncert_resp_mg_C_m3_d)
summary(R.on.D.slr)

errbar(y = Depth.T.BP.Resp.Correlations$Resp_mg_C_m3_d, x = log(Depth.T.BP.Resp.Correlations$Depth_Resp), col=col.Resp, pch=pch.Resp, bg=bg.Resp, ylab="Water col. respiration rate", xlab="log (Depth)", yplus = Depth.T.BP.Resp.Correlations$Resp_mg_C_m3_d+Depth.T.BP.Resp.Correlations$Uncert_resp_mg_C_m3_d, yminus = Depth.T.BP.Resp.Correlations$Resp_mg_C_m3_d-Depth.T.BP.Resp.Correlations$Uncert_resp_mg_C_m3_d)
abline(R.on.D.slr)

pdf(file = "R.on.D.pdf",
    width = 2, height = 2, pointsize = 9,
    bg = "white")
par(pty = "s", mar=c(4,4,1,1))

errbar(y = Depth.T.BP.Resp.Correlations$Resp_mg_C_m3_d, x = log(Depth.T.BP.Resp.Correlations$Depth_Resp), col=col.Resp, pch=pch.Resp, bg=bg.Resp, ylab="Water col. respiration rate", xlab="log (Depth)", yplus = Depth.T.BP.Resp.Correlations$Resp_mg_C_m3_d+Depth.T.BP.Resp.Correlations$Uncert_resp_mg_C_m3_d, yminus = Depth.T.BP.Resp.Correlations$Resp_mg_C_m3_d-Depth.T.BP.Resp.Correlations$Uncert_resp_mg_C_m3_d,xaxt="n",yaxt="n")

axis(1, tck=-0.05)
axis(2, tck=-0.05) 
abline(R.on.D.slr)
dev.off()

# specify symbols and colors for BP
Sta_BP <- factor(Depth.T.BP.Resp.Correlations$Sta_BP)
pch.BP <- c(rep(21,6),rep(22,6),rep(14,6),rep(24,6),rep(25,6),rep(23,6))
col.BP <- c(1)
bg.BP <- c(rep(rgb(80,108,179,maxColorValue=255),6),
           rep(rgb(194,206,233,maxColorValue=255),6),
           rep(rgb(100,14,14,maxColorValue=255),6),
           rep(rgb(153,27,30,maxColorValue=255),6),
           rep(rgb(237,34,36,maxColorValue=255),6),
           rep(rgb(249,191,191,maxColorValue=255),6))

BP.on.T.slr <- lm(Depth.T.BP.Resp.Correlations$BP_mg_C_m3_d ~ Depth.T.BP.Resp.Correlations$Temp_BP, weights = Depth.T.BP.Resp.Correlations$Uncert_BP_mg_C_m3_d)
summary(BP.on.T.slr)

errbar(y = Depth.T.BP.Resp.Correlations$BP_mg_C_m3_d, x = Depth.T.BP.Resp.Correlations$Temp_BP, col=col.BP, pch=pch.BP, bg=bg.BP, ylab="Water col. bacterial production", xlab="Temperature", yplus = Depth.T.BP.Resp.Correlations$BP_mg_C_m3_d+ Depth.T.BP.Resp.Correlations$Uncert_BP_mg_C_m3_d, yminus = Depth.T.BP.Resp.Correlations$BP_mg_C_m3_d-Depth.T.BP.Resp.Correlations$Uncert_BP_mg_C_m3_d)
abline(BP.on.T.slr)

BP.on.D.slr <- lm(Depth.T.BP.Resp.Correlations$BP_mg_C_m3_d ~ log(Depth.T.BP.Resp.Correlations$Depth_BP), weights = Depth.T.BP.Resp.Correlations$Uncert_BP_mg_C_m3_d)
summary(BP.on.D.slr)

errbar(y = Depth.T.BP.Resp.Correlations$BP_mg_C_m3_d, x = log(Depth.T.BP.Resp.Correlations$Depth_BP), col=col.BP, pch=pch.BP, bg=bg.BP, ylab="Water col. bacterial production", xlab="log (Depth)", yplus = Depth.T.BP.Resp.Correlations$BP_mg_C_m3_d+Depth.T.BP.Resp.Correlations$Uncert_BP_mg_C_m3_d, yminus = Depth.T.BP.Resp.Correlations$BP_mg_C_m3_d-Depth.T.BP.Resp.Correlations$Uncert_BP_mg_C_m3_d)
abline(BP.on.D.slr)

pdf(file = "BP.on.D.pdf",
    width = 2, height = 2, pointsize = 9,
    bg = "white")
par(pty = "s", mar=c(4,4,1,1))
errbar(y = Depth.T.BP.Resp.Correlations$BP_mg_C_m3_d, x = log(Depth.T.BP.Resp.Correlations$Depth_BP), col=col.BP, pch=pch.BP, bg=bg.BP, ylab="Water col. bacterial production", xlab="log (Depth)", yplus = Depth.T.BP.Resp.Correlations$BP_mg_C_m3_d+Depth.T.BP.Resp.Correlations$Uncert_BP_mg_C_m3_d, yminus = Depth.T.BP.Resp.Correlations$BP_mg_C_m3_d-Depth.T.BP.Resp.Correlations$Uncert_BP_mg_C_m3_d,xaxt="n",yaxt="n")
axis(1, tck=-0.05)
axis(2, tck=-0.05) 
abline(BP.on.D.slr)
dev.off()

# specify symbols and colors for BP on R
Sta_Comp <- factor(Depth.T.BP.Resp.Correlations$Sta_Comp)
pch.Comp <- c(rep(25,4),rep(23,4),rep(24,5),21,22,14,23)
col.Comp <- c(1)
bg.Comp <- c(rep(rgb(237,34,36,maxColorValue=255),4),
             rep(rgb(249,191,191,maxColorValue=255),4),
             rep(rgb(153,27,30,maxColorValue=255),5),
             rgb(80,108,179,maxColorValue=255),
             rgb(194,206,233,maxColorValue=255),
			 rgb(100,14,14,maxColorValue=255),
             rgb(249,191,191,maxColorValue=255)           
             )

# ***** excluding from the regression one observation (11; 1.627410,0.1385) which is an obvious outlier

BP.on.R.slr <- lm(log(Depth.T.BP.Resp.Correlations$BP_mg_C_m3_d_comp[-11])~log(Depth.T.BP.Resp.Correlations$Resp_mg_C_m3_d_comp[-11]), weights = log(Depth.T.BP.Resp.Correlations$Uncert_resp_mg_C_m3_d_comp[-11]))
summary(BP.on.R.slr)
plot(log(Depth.T.BP.Resp.Correlations$BP_mg_C_m3_d_comp)~log(Depth.T.BP.Resp.Correlations$Resp_mg_C_m3_d_comp), col=col.Comp, pch=pch.Comp, bg=bg.Comp, xlab="log (Water col. respiration rate)", ylab="log (Water col. BP)",xaxt="n",yaxt="n")
abline(BP.on.R.slr)

pdf(file = "BP.on.R.slr.pdf",
    width = 2, height = 2, pointsize = 9,
    bg = "white")
par(pty = "s", mar=c(4,4,1,1))
plot(log(Depth.T.BP.Resp.Correlations$BP_mg_C_m3_d_comp)~log(Depth.T.BP.Resp.Correlations$Resp_mg_C_m3_d_comp), col=col.Comp, pch=pch.Comp, bg=bg.Comp, xlab="log (Water col. respiration rate)", ylab="log (Water col. BP)",xaxt="n",yaxt="n")
axis(1, tck=-0.05)
axis(2, tck=-0.05) 
abline(BP.on.R.slr)

dev.off()


BP.on.R.lrmod2 <- lmodel2(log(Depth.T.BP.Resp.Correlations$BP_mg_C_m3_d_comp[-11]) ~ log(Depth.T.BP.Resp.Correlations$Resp_mg_C_m3_d_comp[-11]),nperm = 99)
plot(BP.on.R.lrmod2, col=col.Comp, pch=pch.Comp, bg=bg.Comp, ylab="log (Water col. bacterial production)", xlab="log (Water col. respiration rate)",xaxt="n",yaxt="n",xlim=c(0.5,4.1))
points(log(Depth.T.BP.Resp.Correlations$Resp_mg_C_m3_d_comp),log(Depth.T.BP.Resp.Correlations$BP_mg_C_m3_d_comp), col=col.Comp, pch=pch.Comp, bg=bg.Comp)
BP.on.R.lrmod2

pdf(file = "BP.on.R.lrmod2.pdf",
    width = 2, height = 2, pointsize = 9,
    bg = "white")
par(pty = "s", mar=c(4,4,1,1))
plot(BP.on.R.lrmod2, col=col.Comp, pch=pch.Comp, bg=bg.Comp, ylab="log (Water col. bacterial production)", xlab="log (Water col. respiration rate)",xaxt="n",yaxt="n",xlim=c(0.1,4.1))
points(log(Depth.T.BP.Resp.Correlations$Resp_mg_C_m3_d_comp),log(Depth.T.BP.Resp.Correlations$BP_mg_C_m3_d_comp), col=col.Comp, pch=pch.Comp, bg=bg.Comp)
axis(1, tck=-0.05)
axis(2, tck=-0.05) 
dev.off()

