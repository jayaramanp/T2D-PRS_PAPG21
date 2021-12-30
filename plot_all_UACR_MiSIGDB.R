#!/bin/R

library(dplyr)
library(ggplot2)

pdf("UACR_CP_densityplot.pdf", width=10, height=10)

print("CP")
file="UACR_DM_PRSet_PAPG21-3_1000gSAS_GRCh37.c2CPv74.best"
group="KEGG_AXON_GUIDANCE"

magic_CP = read.table(file, header=T, sep=" ", stringsAsFactor=F)
myval_CP = magic_CP[,group][468]
sample_PRS.mean = mean(magic_CP[,group])
#tiff("MAGIC_CP_densityplot.tiff", units="in", width=10, height=10, res=300)
ggplot(magic_CP, aes_string(x=group)) + geom_density() + geom_vline(data=magic_CP, aes(xintercept=sample_PRS.mean), linetype="dashed", size=1) + geom_vline(data=magic_CP, aes(xintercept=myval_CP), color = "blue", size = 0.5) + 
labs(x="PRS SCORES", y="DENSITY", title = group)

print("CGP")
##LIU_IL13_MEMORY_MODEL_UP for MAGIC_FG_PRSet_PAPG21-3_1000gSAS_GRCh37.c2CGPv74.best
file="UACR_DM_PRSet_PAPG21-3_1000gSAS_GRCh37.c2CGPv74.best"
group="SMID_BREAST_CANCER_RELAPSE_IN_BONE_UP"
magic_CP = read.table(file, header=T, sep=" ", stringsAsFactor=F)
myval_CP = magic_CP[,group][468]
sample_PRS.mean = mean(magic_CP[,group])
#tiff("MAGIC_CP_densityplot.tiff", units="in", width=10, height=10, res=300)
ggplot(magic_CP, aes_string(x=group)) + geom_density() + geom_vline(data=magic_CP, aes(xintercept=sample_PRS.mean), linetype="dashed", size=1) + geom_vline(data=magic_CP, aes(xintercept=myval_CP), color = "blue", size = 0.5) +
labs(x="PRS SCORES", y="DENSITY", title = group)

print("hALL")
#HALLMARK_UNFOLDED_PROTEIN_RESPONSE for 
file="UACR_DM_PRSet_PAPG21-3_1000gSAS_GRCh37.hALLv74.best"
group="HALLMARK_UNFOLDED_PROTEIN_RESPONSE"
magic_CP = read.table(file, header=T, sep=" ", stringsAsFactor=F)
myval_CP = magic_CP[,group][468]
sample_PRS.mean = mean(magic_CP[,group])
#tiff("MAGIC_CP_densityplot.tiff", units="in", width=10, height=10, res=300)
ggplot(magic_CP, aes_string(x=group)) + geom_density() + geom_vline(data=magic_CP, aes(xintercept=sample_PRS.mean), linetype="dashed", size=1) + geom_vline(data=magic_CP, aes(xintercept=myval_CP), color = "blue", size = 0.5) +
labs(x="PRS SCORES", y="DENSITY", title = group)

print("REGALL")
file="UACR_DM_PRSet_PAPG21-3_1000gSAS_GRCh37.c3REGALLv74.best"
group="AP1_Q2_01"
magic_CP = read.table(file, header=T, sep=" ", stringsAsFactor=F)
myval_CP = magic_CP[,group][468]
sample_PRS.mean = mean(magic_CP[,group])
#tiff("MAGIC_CP_densityplot.tiff", units="in", width=10, height=10, res=300)
ggplot(magic_CP, aes_string(x=group)) + geom_density() + geom_vline(data=magic_CP, aes(xintercept=sample_PRS.mean), linetype="dashed", size=1) + geom_vline(data=magic_CP, aes(xintercept=myval_CP), color = "blue", size = 0.5) +
labs(x="PRS SCORES", y="DENSITY", title = group)


print("ALL")
file="UACR_DM_PRSet_PAPG21-3_1000gSAS_GRCh37.c5ALLv74.best"
group="GOMF_OXIDOREDUCTASE_ACTIVITY_ACTING_ON_NAD_P_H_OXYGEN_AS_ACCEPTOR"
magic_CP = read.table(file, header=T, sep=" ", stringsAsFactor=F)
myval_CP = magic_CP[,group][468]
sample_PRS.mean = mean(magic_CP[,group])
#tiff("MAGIC_CP_densityplot.tiff", units="in", width=10, height=10, res=300)
ggplot(magic_CP, aes_string(x=group)) + geom_density() + geom_vline(data=magic_CP, aes(xintercept=sample_PRS.mean), linetype="dashed", size=1) + geom_vline(data=magic_CP, aes(xintercept=myval_CP), color = "blue", size = 0.5) +
labs(x="PRS SCORES", y="DENSITY", title = group)


dev.off()



