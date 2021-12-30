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

#--base METAANALYSIS_DIAGRAM_SE1.BMI.MOD.nodup.QC.txt --extract DIAGRAM2017_META_PRSice_avg_PAPG21-3_1000gSAS_GRCh37.valid --target PAPG21-3_capstone_1000g_allphase3SAS_GRCh37.QC-merged_excludemissnp.final --no-default --chr CHR --A1 Allele1 --A2 Allele2 --snp SNP --chr-id c:l --pvalue P-value --bp BP --beta --cov PAPG21-3_capstone1000G_merged.QC.COV.final  --pheno PAPG21-3_capstone_1000g_allphase3SAS_GRCh37.QC-merged_excludemissnp.pheno.final --pheno-col T2D --ignore-fid --missing CENTER --thread 2  --stat Effect

Rscript PRSice.R --dir . --prsice PRSice_linux --base METAANALYSIS_DIAGRAM_SE1.BMI.MOD.nodup.QC.txt --extract DIAGRAM2017_META_PRSice_avg_PAPG21-3_1000gSAS_GRCh37.valid --target PAPG21-3_capstone_1000g_allphase3SAS_GRCh37.QC-merged_excludemissnp.final --no-default --chr CHR --A1 Allele1 --A2 Allele2 --snp SNP --chr-id c:l --pvalue P-value --bp BP --beta --cov PAPG21-3_capstone1000G_merged.QC.COV.final  --pheno PAPG21-3_capstone_1000g_allphase3SAS_GRCh37.QC-merged_excludemissnp.pheno.final --pheno-col T2D --ignore-fid --missing CENTER --thread 2  --stat Effect  --gtf misigD_genesets/Homo_sapiens.GRCh37.87.gtf --binary-target T  --msigdb misigD_genesets/c2.cp.v7.4.symbols.gmt --multi-plot 10  --out DIAGRAM_SE1_PRSet_PAPG21-3_1000gSAS_GRCh37.c2CPv74

Rscript PRSice.R --dir . --prsice PRSice_linux --base METAANALYSIS_DIAGRAM_SE1.BMI.MOD.nodup.QC.txt --extract DIAGRAM2017_META_PRSice_avg_PAPG21-3_1000gSAS_GRCh37.valid --target PAPG21-3_capstone_1000g_allphase3SAS_GRCh37.QC-merged_excludemissnp.final --no-default --chr CHR --A1 Allele1 --A2 Allele2 --snp SNP --chr-id c:l --pvalue P-value --bp BP --beta --cov PAPG21-3_capstone1000G_merged.QC.COV.final  --pheno PAPG21-3_capstone_1000g_allphase3SAS_GRCh37.QC-merged_excludemissnp.pheno.final --pheno-col T2D --ignore-fid --missing CENTER --thread 2  --stat Effect  --gtf misigD_genesets/Homo_sapiens.GRCh37.87.gtf --binary-target T  --msigdb misigD_genesets/c2.cgp.v7.4.symbols.gmt --multi-plot 10  --out DIAGRAM_SE1_PRSet_PAPG21-3_1000gSAS_GRCh37.c2CGPv74

Rscript PRSice.R --dir . --prsice PRSice_linux --base METAANALYSIS_DIAGRAM_SE1.BMI.MOD.nodup.QC.txt --extract DIAGRAM2017_META_PRSice_avg_PAPG21-3_1000gSAS_GRCh37.valid --target PAPG21-3_capstone_1000g_allphase3SAS_GRCh37.QC-merged_excludemissnp.final --no-default --chr CHR --A1 Allele1 --A2 Allele2 --snp SNP --chr-id c:l --pvalue P-value --bp BP --beta --cov PAPG21-3_capstone1000G_merged.QC.COV.final  --pheno PAPG21-3_capstone_1000g_allphase3SAS_GRCh37.QC-merged_excludemissnp.pheno.final --pheno-col T2D --ignore-fid --missing CENTER --thread 2  --stat Effect  --gtf misigD_genesets/Homo_sapiens.GRCh37.87.gtf --binary-target T  --msigdb misigD_genesets/h.all.v7.4.symbols.gmt --multi-plot 10  --out DIAGRAM_SE1_PRSet_PAPG21-3_1000gSAS_GRCh37.hALLv74

Rscript PRSice.R --dir . --prsice PRSice_linux --base METAANALYSIS_DIAGRAM_SE1.BMI.MOD.nodup.QC.txt --extract DIAGRAM2017_META_PRSice_avg_PAPG21-3_1000gSAS_GRCh37.valid --target PAPG21-3_capstone_1000g_allphase3SAS_GRCh37.QC-merged_excludemissnp.final --no-default --chr CHR --A1 Allele1 --A2 Allele2 --snp SNP --chr-id c:l --pvalue P-value --bp BP --beta --cov PAPG21-3_capstone1000G_merged.QC.COV.final  --pheno PAPG21-3_capstone_1000g_allphase3SAS_GRCh37.QC-merged_excludemissnp.pheno.final --pheno-col T2D --ignore-fid --missing CENTER --thread 2  --stat Effect  --gtf misigD_genesets/Homo_sapiens.GRCh37.87.gtf --binary-target T  --msigdb misigD_genesets/c3.regall.v7.4.symbols.gmt --multi-plot 10  --out DIAGRAM_SE1_PRSet_PAPG21-3_1000gSAS_GRCh37.c3REGALLv74

Rscript PRSice.R --dir . --prsice PRSice_linux --base METAANALYSIS_DIAGRAM_SE1.BMI.MOD.nodup.QC.txt --extract DIAGRAM2017_META_PRSice_avg_PAPG21-3_1000gSAS_GRCh37.valid --target PAPG21-3_capstone_1000g_allphase3SAS_GRCh37.QC-merged_excludemissnp.final --no-default --chr CHR --A1 Allele1 --A2 Allele2 --snp SNP --chr-id c:l --pvalue P-value --bp BP --beta --cov PAPG21-3_capstone1000G_merged.QC.COV.final  --pheno PAPG21-3_capstone_1000g_allphase3SAS_GRCh37.QC-merged_excludemissnp.pheno.final --pheno-col T2D --ignore-fid --missing CENTER --thread 2  --stat Effect  --gtf misigD_genesets/Homo_sapiens.GRCh37.87.gtf --binary-target T  --msigdb misigD_genesets/c5.all.v7.4.symbols.gmt --multi-plot 10  --out DIAGRAM_SE1_PRSet_PAPG21-3_1000gSAS_GRCh37.c5ALLv74
