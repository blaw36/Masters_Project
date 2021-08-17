library(data.table)

hmt_pvals = readRDS("hmt_pvals.RDS")
nohmt_pvals = readRDS("nohmt_pvals.RDS")

hmt_pvals = hmt_pvals[!is.na(pval)]
nohmt_pvals = nohmt_pvals[!is.na(pval)]

pval = data.table(
  pval.hmt = as.numeric(hmt_pvals$pval)
  ,pval.nohmt = as.numeric(nohmt_pvals$pval)
)

# check histograms
pdf("hmt_hist.pdf")
hist(pval$pval.hmt, breaks=100, main = "hmt")
dev.off()

pdf("nohmt_hist.pdf")
hist(pval$pval.nohmt, breaks=100, main = "nohmt")
dev.off()

# check
# apply qvalue package
library("qvalue")
qval.hmt = qvalue(pval$pval.hmt) # 1
qval.nohmt = qvalue(pval$pval.nohmt) # 0.002288085

# check the proportion of null cases
qval.hmt$pi0
qval.nohmt$pi0

# possible values of FDR
alpha.list = seq(0, 0.1, by=0.001)
length(alpha.list)
## 101

# count the number of significant tests at a given FDR
num.hmt = num.nohmt = rep(NA, length(alpha.list))
for(i in 1:length(alpha.list)){
  num.hmt[i] = sum(qval.hmt$qvalues < alpha.list[i])
  num.nohmt[i] = sum(qval.nohmt$qvalues < alpha.list[i])
}

# number of significant tests at FDR = 0.05
wh = which(alpha.list == 0.05)
num.hmt[wh] # 43
num.nohmt[wh] # 9410 (all of them!)

# Make FDR curves
pdf("fdr_curve.pdf")
hmt.col = "#483D8B"
nohmt.col = "#FF8C00"
par(mar = c(4, 4, 1, 1))
ymax = max(num.hmt, num.nohmt) + 50
ymin = 0
plot(alpha.list, num.hmt, ylim=c(ymin,ymax), col=hmt.col, type = "l", lty = 1, lwd = 1.5, xlab = "FDR", ylab="number of significant tests", main="")
points(alpha.list, num.nohmt, ylim=c(ymin,ymax), col=nohmt.col, type="l", lty = 1, lwd = 1.5)
abline(v=0.05, col="grey")
legend(0,ymax, c("hmt", "nohmt"), col = c(hmt.col, nohmt.col), lty = c(1,1), cex = 0.9, lwd = 1.5, text.col = "black",merge = FALSE, bg = "white")
dev.off()

# if you want to make smooth lines...
wh = max(which(num.hmt == 0))
num.hmt[1:wh] = num.hmt[wh+1]*seq(0,wh-1)/wh 
wh = max(which(num.nohmt == 0))
num.nohmt[1:wh] = num.nohmt[wh+1]*seq(0,wh-1)/wh 

pdf("fdr_curve_smooth.pdf")
hmt.col = "#483D8B"
nohmt.col = "#FF8C00"
par(mar = c(4, 4, 1, 1))
ymax = max(num.hmt, num.nohmt) + 50
ymin = 0
plot(alpha.list, num.hmt, ylim=c(ymin,ymax), col=hmt.col, type = "l", lty = 1, lwd = 1.5, xlab = "FDR", ylab="number of significant tests", main="")
points(alpha.list, num.nohmt, ylim=c(ymin,ymax), col=nohmt.col, type="l", lty = 1, lwd = 1.5)
abline(v=0.05, col="grey")
legend(0,ymax, c("hmt", "nohmt"), col = c(hmt.col, nohmt.col), lty = c(1,1), cex = 0.9, lwd = 1.5, text.col = "black",merge = FALSE, bg = "white")
dev.off()
