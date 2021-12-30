#!/bin/R

library(dplyr)
library(ggplot2)

pdf("MAGIC_CP_densityplot.pdf", width=10, height=10)

print("CP")
#REACTOME_E3_UBIQUITIN_LIGASES_UBIQUITINATE_TARGET_PROTEINS
magic_CP = read.table("MAGIC_FG_PRSet_PAPG21-3_1000gSAS_GRCh37.c2CPv74.best", header=T, sep=" ", stringsAsFactor=F)
myval_CP = magic_CP$REACTOME_E3_UBIQUITIN_LIGASES_UBIQUITINATE_TARGET_PROTEINS[468]
sample_PRS.mean = mean(magic_CP$REACTOME_E3_UBIQUITIN_LIGASES_UBIQUITINATE_TARGET_PROTEINS)
#tiff("MAGIC_CP_densityplot.tiff", units="in", width=10, height=10, res=300)
ggplot(magic_CP, aes(x=REACTOME_E3_UBIQUITIN_LIGASES_UBIQUITINATE_TARGET_PROTEINS)) + geom_density() + geom_vline(data=magic_CP, aes(xintercept=sample_PRS.mean), linetype="dashed", size=1) + geom_vline(data=magic_CP, aes(xintercept=myval_CP), color = "blue", size = 0.5) + 
labs(x="PRS SCORES", y="DENSITY", title = "REACTOME_E3_UBIQUITIN_LIGASES_UBIQUITINATE_TARGET_PROTEINS")

print("CGP")
group = "TCGA_GLIOBLASTOMA_COPY_NUMBER_DN"
##LIU_IL13_MEMORY_MODEL_UP for MAGIC_FG_PRSet_PAPG21-3_1000gSAS_GRCh37.c2CGPv74.best
magic_CGP = read.table("MAGIC_FG_PRSet_PAPG21-3_1000gSAS_GRCh37.c2CGPv74.best", header=T, sep=" ", stringsAsFactor=F)
myval_CGP = magic_CGP[, group][468]
sample_PRS.mean = mean(magic_CGP[,group])
#tiff("MAGIC_CP_densityplot.tiff", units="in", width=10, height=10, res=300)
ggplot(magic_CGP, aes_string(x=group)) + geom_density() + geom_vline(data=magic_CGP, aes(xintercept=sample_PRS.mean), linetype="dashed", size=1) + geom_vline(data=magic_CGP, aes(xintercept=myval_CGP), color = "blue", size = 0.5) +
labs(x="PRS SCORES", y="DENSITY", title = group)

print("hALL")
#HALLMARK_UNFOLDED_PROTEIN_RESPONSE for 
file="MAGIC_FG_PRSet_PAPG21-3_1000gSAS_GRCh37.hALLv74.best"
group="HALLMARK_ALLOGRAFT_REJECTION"
magic_hall = read.table(file, header=T, sep=" ", stringsAsFactor=F)
myval_hall = dplyr::last(magic_hall[,group])
print(myval_hall)
sample_PRS.mean = mean(magic_hall[,group])
sample_PRS.mean

#tiff("MAGIC_CP_densityplot.tiff", units="in", width=10, height=10, res=300)
plot_hall = ggplot(magic_hall, aes_string(x=group)) 
plot_hall + geom_density() + 
geom_vline(data=magic_hall, aes(xintercept=sample_PRS.mean), linetype="dashed", size=1) + 
geom_vline(data=magic_hall, aes(xintercept=myval_hall), color = "blue", size = 0.5) + 
labs(x="PRS SCORES", y="DENSITY", title=group)


print("REGALL")
file="MAGIC_FG_PRSet_PAPG21-3_1000gSAS_GRCh37.c3REGALLv74.best"
group="MIR4310"
magic_regall = read.table(file, header=T, sep=" ", stringsAsFactor=F)
myval_regall = magic_regall[,group][468]
sample_PRS.mean = mean(magic_regall[,group])
#tiff("MAGIC_CP_densityplot.tiff", units="in", width=10, height=10, res=300)
ggplot(magic_regall, aes_string(x=group)) + geom_density() + geom_vline(data=magic_regall, aes(xintercept=sample_PRS.mean), linetype="dashed", size=1) + geom_vline(data=magic_regall, aes(xintercept=myval_regall), color = "blue", size = 0.5) +
labs(x="PRS SCORES", y="DENSITY", title = group)


print("ALL")
file="MAGIC_FG_PRSet_PAPG21-3_1000gSAS_GRCh37.c5ALLv74.best"
group="HP_NEOPLASM_OF_THE_MALE_EXTERNAL_GENITALIA"
magic_all = read.table(file, header=T, sep=" ", stringsAsFactor=F)
myval_all = magic_all[,group][468]
sample_PRS.mean = mean(magic_all[,group])
#tiff("MAGIC_CP_densityplot.tiff", units="in", width=10, height=10, res=300)
ggplot(magic_all, aes_string(x=group)) + geom_density() + geom_vline(data=magic_all, aes(xintercept=sample_PRS.mean), linetype="dashed", size=1) + geom_vline(data=magic_all, aes(xintercept=myval_all), color = "blue", size = 0.5) +
labs(x="PRS SCORES", y="DENSITY", title = group)


dev.off()



