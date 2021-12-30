#!/bin/R

library(dplyr)
library(ggplot2)

pdf("ALL_MAGICDIAGRAMT2DUKBBUACR_PRS_densityplot.pdf", width=10, height=10)
phenofile="PAPG21-3_capstone_1000g_allphase3SAS_GRCh37.QC-merged_excludemissnp.pheno.final"
pheno = read.table(phenofile, header=T, sep=" ",stringsAsFactor=F)

print("UKBB")
file="UKBB_EXOME_CAD_PRSice_avg_PAPG21-3_1000gSAS_GRCh37.best"
group="PAPG21-3"
magic_all = read.table(file, header=T, sep=" ", stringsAsFactor=F)
magic_all$T2D = pheno$T2D[pheno$IID==magic_all$IID]
magic_all$Pheno = magic_all$T2D
myval_all = magic_all$PRS[468]
sample_PRS.mean = mean(magic_all$PRS)
magic_all = magic_all %>%
  mutate(Pheno = ifelse(as.character(T2D) == "1", "Ctrl", "Case"))
#tiff("MAGIC_CP_densityplot.tiff", units="in", width=10, height=10, res=300)
ggplot(magic_all, aes_string(x="PRS", fill="Pheno")) + geom_density(alpha=0.3, size=1) +
geom_vline(data=magic_all, aes(xintercept=sample_PRS.mean), linetype="dashed", size=1) + geom_vline(data=magic_all, aes(xintercept=myval_all), color = "blue", size = 0.5) +
labs(x="PRS SCORES", y="DENSITY", title = "UKBB")

print("UACR")
file="UACR_DM_20170020_PRSice_avg_PAPG21-3_1000gSAS_GRCh37.best"
group="PAPG21-3"
magic_all = read.table(file, header=T, sep=" ", stringsAsFactor=F)
magic_all$T2D = pheno$T2D[pheno$IID==magic_all$IID]
magic_all$Pheno = magic_all$T2D
myval_all = magic_all$PRS[468]
sample_PRS.mean = mean(magic_all$PRS)
magic_all = magic_all %>%
  mutate(Pheno = ifelse(as.character(T2D) == "1", "Ctrl", "Case"))
#tiff("MAGIC_CP_densityplot.tiff", units="in", width=10, height=10, res=300)
ggplot(magic_all, aes_string(x="PRS", fill="Pheno")) + geom_density(alpha=0.3, size=1) +
geom_vline(data=magic_all, aes(xintercept=sample_PRS.mean), linetype="dashed", size=1) + geom_vline(data=magic_all, aes(xintercept=myval_all), color = "blue", size = 0.5) +
labs(x="PRS SCORES", y="DENSITY", title = "UACR")

print("MAGIC")
file="MAGIC_FG_PRSice_avg_PAPG21-3_1000gSAS_GRCh37.best"
group="PAPG21-3"
magic_all = read.table(file, header=T, sep=" ", stringsAsFactor=F)
magic_all$T2D = pheno$T2D[pheno$IID==magic_all$IID]
magic_all$Pheno = magic_all$T2D
myval_all = magic_all$PRS[468]
sample_PRS.mean = mean(magic_all$PRS)
magic_all = magic_all %>%
  mutate(Pheno = ifelse(as.character(T2D) == "1", "Ctrl", "Case"))
#tiff("MAGIC_CP_densityplot.tiff", units="in", width=10, height=10, res=300)
ggplot(magic_all, aes_string(x="PRS", fill="Pheno")) + geom_density(alpha=0.3, size=1) +
geom_vline(data=magic_all, aes(xintercept=sample_PRS.mean), linetype="dashed", size=1) + geom_vline(data=magic_all, aes(xintercept=myval_all), color = "blue", size = 0.5) +
labs(x="PRS SCORES", y="DENSITY", title = "MAGIC")


print("T2D")
file="T2D_TE_PRSice_avg_PAPG21-3_1000gSAS_GRCh37.best"
group="PAPG21-3"
magic_all = read.table(file, header=T, sep=" ", stringsAsFactor=F)
#head(pheno)
magic_all$T2D = pheno$T2D[pheno$IID==magic_all$IID]
magic_all$Pheno = magic_all$T2D
#magic_all$Pheno = magic_all$T2D[magic_all$T2D == 2] <- "Case"
magic_all = magic_all %>% 
  mutate(Pheno = ifelse(as.character(T2D) == "1", "Ctrl", "Case"))
myval_all = magic_all$PRS[468]
sample_PRS.mean = mean(magic_all$PRS)
#tiff("MAGIC_CP_densityplot.tiff", units="in", width=10, height=10, res=300)
ggplot(magic_all, aes_string(x="PRS", fill="Pheno")) + geom_density(alpha=0.3, size=1) +  
 geom_vline(data=magic_all, aes(xintercept=sample_PRS.mean), linetype="dashed", size=1) + geom_vline(data=magic_all, aes(xintercept=myval_all), color = "blue", size = 0.5) +
labs(x="PRS SCORES", y="DENSITY", title = "T2D")


print("DIAGRAM")
file="DIAGRAM2017_META_PRSice_avg_PAPG21-3_1000gSAS_GRCh37.best"
group="PAPG21-3"
magic_all = read.table(file, header=T, sep=" ", stringsAsFactor=F)
magic_all$T2D = pheno$T2D[pheno$IID==magic_all$IID]
magic_all$Pheno = magic_all$T2D
myval_all = magic_all$PRS[468]
sample_PRS.mean = mean(magic_all$PRS)
magic_all = magic_all %>%
  mutate(Pheno = ifelse(as.character(T2D) == "1", "Ctrl", "Case"))
#tiff("MAGIC_CP_densityplot.tiff", units="in", width=10, height=10, res=300)
ggplot(magic_all, aes_string(x="PRS", fill="Pheno")) + geom_density(alpha=0.3, size=1) +
 geom_vline(data=magic_all, aes(xintercept=sample_PRS.mean), linetype="dashed", size=1) + geom_vline(data=magic_all, aes(xintercept=myval_all), color = "blue", size = 0.5) +
labs(x="PRS SCORES", y="DENSITY", title = "DIAGRAM")

dev.off()

