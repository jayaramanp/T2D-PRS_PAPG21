#!/bin/R

library(dplyr)
library(ggplot2)

pdf("UKBB_CP_densityplot.pdf", width=10, height=10)

print("CP")
file="UKBB_CAD_PRSet_PAPG21-3_1000gSAS_GRCh37.c2CPv74.best"
group="REACTOME_PEROXISOMAL_PROTEIN_IMPORT"

magic_CP = read.table(file, header=T, sep=" ", stringsAsFactor=F)
myval_CP = magic_CP[,group][468]
sample_PRS.mean = mean(magic_CP[,group])
#tiff("MAGIC_CP_densityplot.tiff", units="in", width=10, height=10, res=300)
ggplot(magic_CP, aes_string(x=group)) + geom_density() + geom_vline(data=magic_CP, aes(xintercept=sample_PRS.mean), linetype="dashed", size=1) + geom_vline(data=magic_CP, aes(xintercept=myval_CP), color = "blue", size = 0.5) + 
labs(x="PRS SCORES", y="DENSITY", title = group)

print("CGP")
##LIU_IL13_MEMORY_MODEL_UP for MAGIC_FG_PRSet_PAPG21-3_1000gSAS_GRCh37.c2CGPv74.best
file="UKBB_CAD_PRSet_PAPG21-3_1000gSAS_GRCh37.c2CGPv74.best"
group="SATO_SILENCED_EPIGENETICALLY_IN_PANCREATIC_CANCER"
magic_CP = read.table(file, header=T, sep=" ", stringsAsFactor=F)
myval_CP = magic_CP[,group][468]
sample_PRS.mean = mean(magic_CP[,group])
#tiff("MAGIC_CP_densityplot.tiff", units="in", width=10, height=10, res=300)
ggplot(magic_CP, aes_string(x=group)) + geom_density() + geom_vline(data=magic_CP, aes(xintercept=sample_PRS.mean), linetype="dashed", size=1) + geom_vline(data=magic_CP, aes(xintercept=myval_CP), color = "blue", size = 0.5) +
labs(x="PRS SCORES", y="DENSITY", title = group)

print("hALL")
#HALLMARK_UNFOLDED_PROTEIN_RESPONSE for 
file="UKBB_CAD_PRSet_PAPG21-3_1000gSAS_GRCh37.hALLv74.best"
group="HALLMARK_XENOBIOTIC_METABOLISM"
magic_CP = read.table(file, header=T, sep=" ", stringsAsFactor=F)
myval_CP = magic_CP[,group][468]
sample_PRS.mean = mean(magic_CP[,group])
#tiff("MAGIC_CP_densityplot.tiff", units="in", width=10, height=10, res=300)
ggplot(magic_CP, aes_string(x=group)) + geom_density() + geom_vline(data=magic_CP, aes(xintercept=sample_PRS.mean), linetype="dashed", size=1) + geom_vline(data=magic_CP, aes(xintercept=myval_CP), color = "blue", size = 0.5) +
labs(x="PRS SCORES", y="DENSITY", title = group)

print("REGALL")
file="UKBB_CAD_PRSet_PAPG21-3_1000gSAS_GRCh37.c3REGALLv74.best"
group="MIR3127_5P"
magic_CP = read.table(file, header=T, sep=" ", stringsAsFactor=F)
myval_CP = magic_CP[,group][468]
sample_PRS.mean = mean(magic_CP[,group])
#tiff("MAGIC_CP_densityplot.tiff", units="in", width=10, height=10, res=300)
ggplot(magic_CP, aes_string(x=group)) + geom_density() + geom_vline(data=magic_CP, aes(xintercept=sample_PRS.mean), linetype="dashed", size=1) + geom_vline(data=magic_CP, aes(xintercept=myval_CP), color = "blue", size = 0.5) +
labs(x="PRS SCORES", y="DENSITY", title = group)


print("ALL")
file="UKBB_CAD_PRSet_PAPG21-3_1000gSAS_GRCh37.c5ALLv74.best"
group="GOBP_SURFACTANT_HOMEOSTASIS"
magic_CP = read.table(file, header=T, sep=" ", stringsAsFactor=F)
myval_CP = magic_CP[,group][468]
sample_PRS.mean = mean(magic_CP[,group])
#tiff("MAGIC_CP_densityplot.tiff", units="in", width=10, height=10, res=300)
ggplot(magic_CP, aes_string(x=group)) + geom_density() + geom_vline(data=magic_CP, aes(xintercept=sample_PRS.mean), linetype="dashed", size=1) + geom_vline(data=magic_CP, aes(xintercept=myval_CP), color = "blue", size = 0.5) +
labs(x="PRS SCORES", y="DENSITY", title = group)


dev.off()



