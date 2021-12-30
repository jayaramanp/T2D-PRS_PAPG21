#!bin/R

prs.all <- read.table("UACR_DM_20170020_PRSice_avg_PAPG21-3_1000gSAS_GRCh37.best", header = TRUE)
prs.all <- as.data.frame(prs.all)
# Find the mean
library(plyr)
library(ggplot2)

print("MEAN PRS:")
cprs <- (PRS.mean=mean(prs.all$PRS))
cprs

tiff("UACR_DM_PRS_densityplot.tiff", units="in", width=10, height=10, res=300)
ggplot(prs.all, aes(x=PRS)) + geom_density() + geom_vline(data=prs.all, aes(xintercept=PRS.mean), linetype="dashed", size=1) + geom_vline(data=prs.all, aes(xintercept=-0.418322903), color = "blue", size = 0.5)
dev.off()
