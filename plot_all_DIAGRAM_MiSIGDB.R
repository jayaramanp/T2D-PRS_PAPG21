#!/bin/R

library(dplyr)
library(ggplot2)

pdf("DIAGRAM_CP_densityplot.pdf", width=10, height=10)

print("CGP")
file="DIAGRAM_SE1_PRSet_PAPG21-3_1000gSAS_GRCh37.c2CGPv74.best"
group="MORI_PLASMA_CELL_UP"

magic_CP = read.table(file, header=T, sep=" ", stringsAsFactor=F)
myval_CP = magic_CP[,group][468]
sample_PRS.mean = mean(magic_CP[,group])
#tiff("MAGIC_CP_densityplot.tiff", units="in", width=10, height=10, res=300)
ggplot(magic_CP, aes_string(x=group)) + geom_density() + geom_vline(data=magic_CP, aes(xintercept=sample_PRS.mean), linetype="dashed", size=1) + geom_vline(data=magic_CP, aes(xintercept=myval_CP), color = "blue", size = 0.5) + 
labs(x="PRS SCORES", y="DENSITY", title = group)

print("CP")
##LIU_IL13_MEMORY_MODEL_UP for MAGIC_FG_PRSet_PAPG21-3_1000gSAS_GRCh37.c2CGPv74.best
file="DIAGRAM_SE1_PRSet_PAPG21-3_1000gSAS_GRCh37.c2CPv74.best"
group="WP_ID_SIGNALING_PATHWAY"
magic_CP = read.table(file, header=T, sep=" ", stringsAsFactor=F)
myval_CP = magic_CP[,group][468]
sample_PRS.mean = mean(magic_CP[,group])
#tiff("MAGIC_CP_densityplot.tiff", units="in", width=10, height=10, res=300)
ggplot(magic_CP, aes_string(x=group)) + geom_density() + geom_vline(data=magic_CP, aes(xintercept=sample_PRS.mean), linetype="dashed", size=1) + geom_vline(data=magic_CP, aes(xintercept=myval_CP), color = "blue", size = 0.5) +
labs(x="PRS SCORES", y="DENSITY", title = group)

print("hALL")
#HALLMARK_UNFOLDED_PROTEIN_RESPONSE for 
file="DIAGRAM_SE1_PRSet_PAPG21-3_1000gSAS_GRCh37.hALLv74.best"
group="HALLMARK_UNFOLDED_PROTEIN_RESPONSE"
magic_CP = read.table(file, header=T, sep=" ", stringsAsFactor=F)
myval_CP = magic_CP[,group][468]
sample_PRS.mean = mean(magic_CP[,group])
#tiff("MAGIC_CP_densityplot.tiff", units="in", width=10, height=10, res=300)
ggplot(magic_CP, aes_string(x=group)) + geom_density() + geom_vline(data=magic_CP, aes(xintercept=sample_PRS.mean), linetype="dashed", size=1) + geom_vline(data=magic_CP, aes(xintercept=myval_CP), color = "blue", size = 0.5) +
labs(x="PRS SCORES", y="DENSITY", title = group)

print("REGALL")
file="DIAGRAM_SE1_PRSet_PAPG21-3_1000gSAS_GRCh37.c3REGALLv74.best"
group="SMN1_SMN2_TARGET_GENES"
magic_CP = read.table(file, header=T, sep=" ", stringsAsFactor=F)
myval_CP = magic_CP[,group][468]
sample_PRS.mean = mean(magic_CP[,group])
#tiff("MAGIC_CP_densityplot.tiff", units="in", width=10, height=10, res=300)
ggplot(magic_CP, aes_string(x=group)) + geom_density() + geom_vline(data=magic_CP, aes(xintercept=sample_PRS.mean), linetype="dashed", size=1) + geom_vline(data=magic_CP, aes(xintercept=myval_CP), color = "blue", size = 0.5) +
labs(x="PRS SCORES", y="DENSITY", title = group)


print("ALL")
file="DIAGRAM_SE1_PRSet_PAPG21-3_1000gSAS_GRCh37.c5ALLv74.best"
group="HP_FOCAL_MYOCLONIC_SEIZURE"
magic_CP = read.table(file, header=T, sep=" ", stringsAsFactor=F)
myval_CP = magic_CP[,group][468]
sample_PRS.mean = mean(magic_CP[,group])
#tiff("MAGIC_CP_densityplot.tiff", units="in", width=10, height=10, res=300)
ggplot(magic_CP, aes_string(x=group)) + geom_density() + geom_vline(data=magic_CP, aes(xintercept=sample_PRS.mean), linetype="dashed", size=1) + geom_vline(data=magic_CP, aes(xintercept=myval_CP), color = "blue", size = 0.5) +
labs(x="PRS SCORES", y="DENSITY", title = group)


dev.off()



