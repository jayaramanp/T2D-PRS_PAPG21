#!/bin/R


module load R
module load plink
module load plink2

#######
#-rw-r--r--  1 jayarp02 bsr2401  48K Dec 27 12:24 h.all.v7.4.symbols.gmt
#-rw-r--r--  1 jayarp02 bsr2401 2.6M Dec 27 12:24 c2.cgp.v7.4.symbols.gmt
#-rw-r--r--  1 jayarp02 bsr2401 1.2M Dec 27 12:24 c2.cp.v7.4.symbols.gmt
#-rw-r--r--  1 jayarp02 bsr2401 5.5M Dec 27 12:24 c3.regall.v7.4.symbols.gmt
#-rw-r--r--  1 jayarp02 bsr2401 8.9M Dec 27 12:24 c5.all.v7.4.symbols.gmt
#######

#--base T2D_TranEthnic.BMIadjusted.nodup.QC.final.txt --extract T2D_TE_PRSice_avg_PAPG21-3_1000gSAS_GRCh37.valid --target PAPG21-3_capstone_1000g_allphase3SAS_GRCh37.QC-merged_excludemissnp.final --no-default --chr CHR --A1 A1 --A2 A2 --snp SNP --chr-id c:l --pvalue P --bp BP --beta --cov PAPG21-3_capstone1000G_merged.QC.COV.final  --pheno PAPG21-3_capstone_1000g_allphase3SAS_GRCh37.QC-merged_excludemissnp.pheno.final --pheno-col T2D --ignore-fid --missing CENTER --thread 2  --stat BETA

Rscript PRSice.R --dir . --prsice PRSice_linux --base T2D_TranEthnic.BMIadjusted.nodup.QC.final.txt --extract T2D_TE_PRSice_avg_PAPG21-3_1000gSAS_GRCh37.valid --target PAPG21-3_capstone_1000g_allphase3SAS_GRCh37.QC-merged_excludemissnp.final --no-default --chr CHR --A1 A1 --A2 A2 --snp SNP --chr-id c:l --pvalue P --bp BP --beta --cov PAPG21-3_capstone1000G_merged.QC.COV.final  --pheno PAPG21-3_capstone_1000g_allphase3SAS_GRCh37.QC-merged_excludemissnp.pheno.final --pheno-col T2D --ignore-fid --missing CENTER --thread 2  --stat BETA  --gtf misigD_genesets/Homo_sapiens.GRCh37.87.gtf --binary-target T  --msigdb misigD_genesets/c2.cp.v7.4.symbols.gmt --multi-plot 10  --out T2D_TranEthnic_PRSet_PAPG21-3_1000gSAS_GRCh37.c2CPv74

Rscript PRSice.R --dir . --prsice PRSice_linux --base T2D_TranEthnic.BMIadjusted.nodup.QC.final.txt --extract T2D_TE_PRSice_avg_PAPG21-3_1000gSAS_GRCh37.valid --target PAPG21-3_capstone_1000g_allphase3SAS_GRCh37.QC-merged_excludemissnp.final --no-default --chr CHR --A1 A1 --A2 A2 --snp SNP --chr-id c:l --pvalue P --bp BP --beta --cov PAPG21-3_capstone1000G_merged.QC.COV.final  --pheno PAPG21-3_capstone_1000g_allphase3SAS_GRCh37.QC-merged_excludemissnp.pheno.final --pheno-col T2D --ignore-fid --missing CENTER --thread 2  --stat BETA  --gtf misigD_genesets/Homo_sapiens.GRCh37.87.gtf --binary-target T  --msigdb misigD_genesets/c2.cgp.v7.4.symbols.gmt --multi-plot 10  --out T2D_TranEthnic_PRSet_PAPG21-3_1000gSAS_GRCh37.c2CGPv74

Rscript PRSice.R --dir . --prsice PRSice_linux --base T2D_TranEthnic.BMIadjusted.nodup.QC.final.txt --extract T2D_TE_PRSice_avg_PAPG21-3_1000gSAS_GRCh37.valid --target PAPG21-3_capstone_1000g_allphase3SAS_GRCh37.QC-merged_excludemissnp.final --no-default --chr CHR --A1 A1 --A2 A2 --snp SNP --chr-id c:l --pvalue P --bp BP --beta --cov PAPG21-3_capstone1000G_merged.QC.COV.final  --pheno PAPG21-3_capstone_1000g_allphase3SAS_GRCh37.QC-merged_excludemissnp.pheno.final --pheno-col T2D --ignore-fid --missing CENTER --thread 2  --stat BETA  --gtf misigD_genesets/Homo_sapiens.GRCh37.87.gtf --binary-target T  --msigdb misigD_genesets/h.all.v7.4.symbols.gmt --multi-plot 10  --out T2D_TranEthnic_PRSet_PAPG21-3_1000gSAS_GRCh37.hALLv74

Rscript PRSice.R --dir . --prsice PRSice_linux --base T2D_TranEthnic.BMIadjusted.nodup.QC.final.txt --extract T2D_TE_PRSice_avg_PAPG21-3_1000gSAS_GRCh37.valid --target PAPG21-3_capstone_1000g_allphase3SAS_GRCh37.QC-merged_excludemissnp.final --no-default --chr CHR --A1 A1 --A2 A2 --snp SNP --chr-id c:l --pvalue P --bp BP --beta --cov PAPG21-3_capstone1000G_merged.QC.COV.final  --pheno PAPG21-3_capstone_1000g_allphase3SAS_GRCh37.QC-merged_excludemissnp.pheno.final --pheno-col T2D --ignore-fid --missing CENTER --thread 2  --stat BETA  --gtf misigD_genesets/Homo_sapiens.GRCh37.87.gtf --binary-target T  --msigdb misigD_genesets/c3.regall.v7.4.symbols.gmt --multi-plot 10  --out T2D_TranEthnic_PRSet_PAPG21-3_1000gSAS_GRCh37.c3REGALLv74

Rscript PRSice.R --dir . --prsice PRSice_linux --base T2D_TranEthnic.BMIadjusted.nodup.QC.final.txt --extract T2D_TE_PRSice_avg_PAPG21-3_1000gSAS_GRCh37.valid --target PAPG21-3_capstone_1000g_allphase3SAS_GRCh37.QC-merged_excludemissnp.final --no-default --chr CHR --A1 A1 --A2 A2 --snp SNP --chr-id c:l --pvalue P --bp BP --beta --cov PAPG21-3_capstone1000G_merged.QC.COV.final  --pheno PAPG21-3_capstone_1000g_allphase3SAS_GRCh37.QC-merged_excludemissnp.pheno.final --pheno-col T2D --ignore-fid --missing CENTER --thread 2  --stat BETA  --gtf misigD_genesets/Homo_sapiens.GRCh37.87.gtf --binary-target T  --msigdb misigD_genesets/c5.all.v7.4.symbols.gmt --multi-plot 10  --out T2D_TranEthnic_PRSet_PAPG21-3_1000gSAS_GRCh37.c5ALLv74
