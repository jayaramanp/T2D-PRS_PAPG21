#!/bin/sh

module load R
module load plink
module load plink2


Rscript PRSice.R --dir . --prsice PRSice_linux --base METAANALYSIS_DIAGRAM_SE1.BMI.MOD.nodup.QC.txt --extract DIAGRAM2017_META_PRSice_avg_PAPG21-3_1000gSAS_GRCh37.valid --target PAPG21-3_capstone_1000g_allphase3SAS_GRCh37.QC-merged_excludemissnp.final --no-default --chr CHR --A1 Allele1 --A2 Allele2 --snp SNP --chr-id c:l --pvalue P-value --bp BP --beta --cov PAPG21-3_capstone1000G_merged.QC.COV.final  --pheno PAPG21-3_capstone_1000g_allphase3SAS_GRCh37.QC-merged_excludemissnp.pheno.final --pheno-col T2D --ignore-fid --missing CENTER --thread 2  --stat Effect   --binary-target T   --score avg   --quantile 100  --quant-break 1,5,10,20,40,60,80,90,95,99,100   --quant-ref 60  --out DIAGRAM2017_META_PRSice_avg_PAPG21-3_1000gSAS_GRCh37
