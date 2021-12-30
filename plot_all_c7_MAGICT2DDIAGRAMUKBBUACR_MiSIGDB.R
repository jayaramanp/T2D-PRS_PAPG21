#!/bin/R

library(dplyr)
library(ggplot2)

pdf("ALL_MAGICDIAGRAMT2DUKBBUACR_C7_densityplot.pdf", width=10, height=10)

print("UKBB c7")
file="UKBB_CAD_PRSet_PAPG21-3_1000gSAS_GRCh37.c7ALLv74.best"
group="GSE5589_UNSTIM_VS_45MIN_LPS_STIM_MACROPHAGE_UP"
magic_all = read.table(file, header=T, sep=" ", stringsAsFactor=F)
myval_all = magic_all[,group][468]
sample_PRS.mean = mean(magic_all[,group])
#tiff("MAGIC_CP_densityplot.tiff", units="in", width=10, height=10, res=300)
ggplot(magic_all, aes_string(x=group)) + geom_density() + geom_vline(data=magic_all, aes(xintercept=sample_PRS.mean), linetype="dashed", size=1) + geom_vline(data=magic_all, aes(xintercept=myval_all), color = "blue", size = 0.5) +
labs(x="PRS SCORES", y="DENSITY", title = group)

#UACR_DM_PRSet_PAPG21-3_1000gSAS_GRCh37.c7ALLv74.best

print("UACR c7")
file="UACR_DM_PRSet_PAPG21-3_1000gSAS_GRCh37.c7ALLv74.best"
group="GSE15767_MED_VS_SCS_MAC_LN_DN"
magic_all = read.table(file, header=T, sep=" ", stringsAsFactor=F)
myval_all = magic_all[,group][468]
sample_PRS.mean = mean(magic_all[,group])
#tiff("MAGIC_CP_densityplot.tiff", units="in", width=10, height=10, res=300)
ggplot(magic_all, aes_string(x=group)) + geom_density() + geom_vline(data=magic_all, aes(xintercept=sample_PRS.mean), linetype="dashed", size=1) + geom_vline(data=magic_all, aes(xintercept=myval_all), color = "blue", size = 0.5) +
labs(x="PRS SCORES", y="DENSITY", title = group)

print("MAGIC c7")
file="MAGIC_FG_PRSet_PAPG21-3_1000gSAS_GRCh37.c7ALLv74.best"
group="GSE40274_FOXP3_VS_FOXP3_AND_EOS_TRANSDUCED_ACTIVATED_CD4_TCELL_DN"
magic_all = read.table(file, header=T, sep=" ", stringsAsFactor=F)
myval_all = magic_all[,group][468]
sample_PRS.mean = mean(magic_all[,group])
#tiff("MAGIC_CP_densityplot.tiff", units="in", width=10, height=10, res=300)
ggplot(magic_all, aes_string(x=group)) + geom_density() + geom_vline(data=magic_all, aes(xintercept=sample_PRS.mean), linetype="dashed", size=1) + geom_vline(data=magic_all, aes(xintercept=myval_all), color = "blue", size = 0.5) +
labs(x="PRS SCORES", y="DENSITY", title = group)


print("T2D c7")
file="T2D_TranEthnic_PRSet_PAPG21-3_1000gSAS_GRCh37.c7ALLv74.best"
group="GSE1432_1H_VS_6H_IFNG_MICROGLIA_UP"
magic_all = read.table(file, header=T, sep=" ", stringsAsFactor=F)
myval_all = magic_all[,group][468]
sample_PRS.mean = mean(magic_all[,group])
#tiff("MAGIC_CP_densityplot.tiff", units="in", width=10, height=10, res=300)
ggplot(magic_all, aes_string(x=group)) + geom_density() + geom_vline(data=magic_all, aes(xintercept=sample_PRS.mean), linetype="dashed", size=1) + geom_vline(data=magic_all, aes(xintercept=myval_all), color = "blue", size = 0.5) +
labs(x="PRS SCORES", y="DENSITY", title = group)

print("DIAGRAM c7")
file="DIAGRAM_SE1_PRSet_PAPG21-3_1000gSAS_GRCh37.c7ALLv74.best"
group="GSE3982_EOSINOPHIL_VS_DC_UP"
magic_all = read.table(file, header=T, sep=" ", stringsAsFactor=F)
myval_all = magic_all[,group][468]
sample_PRS.mean = mean(magic_all[,group])
#tiff("MAGIC_CP_densityplot.tiff", units="in", width=10, height=10, res=300)
ggplot(magic_all, aes_string(x=group)) + geom_density() + geom_vline(data=magic_all, aes(xintercept=sample_PRS.mean), linetype="dashed", size=1) + geom_vline(data=magic_all, aes(xintercept=myval_all), color = "blue", size = 0.5) +
labs(x="PRS SCORES", y="DENSITY", title = group)


dev.off()

