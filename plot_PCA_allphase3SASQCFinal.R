#!/bin/R

sas.pca<- read.table("1000G/all_phase3_SAS.QC.final.pca.eigenvec")
sas.pca <- as.data.frame(sas.pca[3:22])
colnames(sas.pca) <- c("PCA_1", "PCA_2", "PCA_3", "PCA_4", "PCA_5", "PCA_6", "PCA_7", "PCA_8", "PCA_9", "PCA_10", "PCA_11", "PCA_12", "PCA_13", "PCA_14", "PCA_15", "PCA_16", "PCA_17", "PCA_18", "PCA_19", "PCA_20")

library(ggfortify) 
pca_res <- prcomp(sas.pca,  scale = TRUE) 
summary(pca_res)
var_explained <-pca_res$sdev^2/sum(pca_res$sdev^2)
var_explained[1:5]
PCA <- as.data.frame(pca_res$x)

tiff("SAS_QC_pca.tiff", units="in", width=10, height=10, res=300)
ggplot(PCA, aes(x=PC1,y=PC2)) + geom_point(size=4) + theme_bw(base_size=24) + labs(x=paste0("PC1: ",round(var_explained[3]*100,1),"%"), y=paste0("PC2: ",round(var_explained[4]*100,1),"%")) + theme(legend.position="top")
dev.off()

covariate <- read.table("1000G/all_phase3_SAS.QC.NEWCOVARIATES", header=T, sep="\t")
colnames(covariate) <- c("FID","IID","Sex","Age","BMI")
pcs <- read.table("1000G/all_phase3_SAS.QC.final.pca.eigenvec", header=F)
colnames(pcs) <- c("FID","IID", paste0("PC",1:6))
cov <- merge(covariate, pcs, by=c("FID", "IID"))
write.table(cov,"1000G/all_phase3_SAS.QC.final.COV", quote=F, row.names=F)
