#!/bin/R

pca_TargetData <- read.table("PAPG21-3_capstone_1000g_allphase3SAS_GRCh37.QC-merged_excludemissnp.pca.eigenvec")
TD_pca <- as.data.frame(pca_TargetData[3:22])
colnames(TD_pca) <- c("PCA_1", "PCA_2", "PCA_3", "PCA_4", "PCA_5", "PCA_6", "PCA_7", "PCA_8", "PCA_9", "PCA_10", "PCA_11", "PCA_12", "PCA_13", "PCA_14", "PCA_15", "PCA_16", "PCA_17", "PCA_18", "PCA_19", "PCA_20")

library(ggfortify) 
pca_res <- prcomp(TD_pca,  scale = TRUE) 
summary(pca_res)
var_explained <-pca_res$sdev^2/sum(pca_res$sdev^2)
var_explained[1:5]
PCA <- as.data.frame(pca_res$x)

tiff("TargetData_PAPGmerged1000G_pca.tiff", units="in", width=10, height=10, res=300)
ggplot(PCA, aes(x=PC1,y=PC2)) + geom_point(color = "blue", size=4) + theme_bw(base_size=24) + labs(x=paste0("PC1: ",round(var_explained[3]*100,1),"%"), y=paste0("PC2: ",round(var_explained[4]*100,1),"%")) + theme(legend.position="top")
dev.off()

covariate <- read.table("PAPG21-3_capstone_1000g_allphase3SAS_GRCh37.QC-merged_excludemissnp.NEWCOVARIATES", header=T, sep="\t")
colnames(covariate) <- c("FID","IID","Sex","Age","BMI")
pcs <- read.table("PAPG21-3_capstone_1000g_allphase3SAS_GRCh37.QC-merged_excludemissnp.pca.eigenvec", header=F)
colnames(pcs) <- c("FID","IID", paste0("PC",1:20))
cov <- merge(covariate, pcs, by=c("FID", "IID"))
write.table(cov,"PAPG21-3_capstone1000G_merged.QC.COV", quote=F, row.names=F)
