# smaller indiv images ----------------------------------------------------

plot.new()
# Open file
# png("images/ch2_high_signal_wc.png", width = 900, height = 600, pointsize = 20)
pdf("images/talk_signal_lowlvl.PDF", width = 6, height = 3, pointsize = 10)
# Create plot
plot(wlet_69$D[1:512]
     , xlim = c(1,512)
     , type = "h"
     , ylab = "WC"
     , xlab = "Location"
     # , main = "Wavelet coefficients"
     , axes = F)
# plot(wlet_69, scaling = "by.level"
#      , ylab = "Scale Level (Resolution)",xlab = "Location"
#      , main = "Wavelet coefficients"
#      , first.level = 9)
axis(1, at = seq(0,512,128),labels = seq(0,512,128))
axis(2, at = 0, labels = 0, tick = FALSE)
# Close the file
dev.off()


# signal image ------------------------------------------------------------
## Load data from WaveQTL

library(plotly)
library(data.table)
library(reshape2)
library(wavethresh)

data.path = "~/Cpp/WaveQTL/data/dsQTL/"
source("~/Cpp/WaveQTL/R/WaveQTL_preprocess_funcs.R")
source("code/sim3_functions.R")

# Read in phenotype data
pheno.dat = as.matrix(read.table(paste0(data.path,"chr17.10160989.10162012.pheno.dat")))

two_datasets <- pheno.dat[c(21,69),]

wlet_21 <- wd(two_datasets[1,], filter.number = 1, family = "DaubExPhase")
wlet_69 <- wd(two_datasets[2,], filter.number = 1, family = "DaubExPhase")

two_datasets_plot = data.frame(t(two_datasets))
two_datasets_plot = melt(two_datasets_plot)
two_datasets_plot$variable = rep(c(21,69), each = 1024)
two_datasets_plot$loc = rep(1:1024, 2)

## Wavelet 69 - more signals
# Plot data and w/let decomps

plot.new()
# Open file
# png("images/ch2_high_signal_seq.png", width = 900, height = 350, pointsize = 20)
pdf("images/talk_high_signal_small.pdf", width = 12, height = 5, pointsize = 20)
# Create plot
plot(two_datasets_plot[two_datasets_plot$variable == 69, "loc"][506:513]
     ,two_datasets_plot[two_datasets_plot$variable == 69, "value"][506:513]
     ,main = ""
     ,xlab = ""
     ,ylab = "Normalised count"
     ,type = "b"
     ,xaxt = "n"
     ,axes = F
     ,ylim = c(0,3)
     ,xlim = c(506,513)
     ,col = "#c00000")
# axis(1, at = seq(506,513))
axis(2)
# Close the file
dev.off()

# Posterior mean e.g: -----------------------------------------------------

# run 'get_effectSizeInDataSpace.R' first

pdf("../../../../Dropbox/Uni Stuff - Masters/Research Project/Masters_Project_Git/images/talk_wqtl_effect_v2.pdf", width = 10, height=5, pointsize = 15)
par(mar = c(4,4,4,2))
plot(1,1,type="n", xlab = "Base Location", ylab = bquote("Effect size, "~alpha[b])
     ,ylim=c(ymin_beta, ymax_beta),xlim=c(1, 1024),main ="Posterior mean +/-3 posterior SDs - LogL: 73.22, p-val: <0.0001", axes=FALSE)
axis(2)
if(length(col_posi) > 0){
        for(j in 1:length(col_posi)){
                polygon(c(col_posi[j]-0.5, col_posi[j]-0.5, col_posi[j]+0.5, col_posi[j]+0.5), c(ymin_beta-2, ymax_beta+2, ymax_beta+2, ymin_beta-2), col ="pink", border = NA)
        }
}

abline(h = 0, col = "red")
points(xval, beta_dataS, col = "blue", type="l")
points(xval, beta_l, col = "skyblue", type="l")
points(xval, beta_r, col = "skyblue", type="l")
box()
axis(side = 1,at = c(1,seq(256,1024,256))
     ,tick = c(1,seq(256,1024,256))
     ,labels = c(1,seq(256,1024,256)))

dev.off()


