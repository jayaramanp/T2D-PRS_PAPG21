#!/bin/R
library(data.table)

df1 = fread("1000G/all_phase3_SAS.QC.fam.nopheno", data.table = F, header=F, stringsAsFactors = F, fill = T) 
head(df1)

df2 = fread("1000G/1kg_phenotype_mock_SAS_2021-12-25.tsv", data.table = F, header=T, stringsAsFactors = F, fill = T)
head(df2)

df3 = df1[, c("V1", "V2", "V5")] 
head(df3)
df3$V3 <- df2$age[match(df1$V2, df2$sample)]
df3$V4 <- df2$bmi[match(df1$V2, df2$sample)]
colnames(df3) <- c("FID","IID","Sex","Age","BMI")
head(df3)

df1$V6 <- df2$t2d[match(df1$V2, df2$sample)]

head(df1)

df4 = df1[,c("V1","V2","V6")]
colnames(df4) <- c("FID", "IID","T2D")
head(df4)

fwrite(df1,
       file = "1000G/all_phase3_SAS.QC.fam", col.names=F,
       sep = "\t", row.names = F)

fwrite(df3,
       file = "1000G/all_phase3_SAS.QC.NEWCOVARIATES", col.names=T,
       sep = "\t", row.names = F)

fwrite(df4,
       file = "1000G/all_phase3_SAS.QC.pheno", col.names=T,
       sep = "\t", row.names = F)

