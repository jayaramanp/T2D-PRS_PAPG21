#!/bin/sh

module load R
module load plink
module load plink2


Rscript PRSice.R --dir . --prsice PRSice_linux --base UKBB.GWAS1KG.EXOME.CAD.SOFT.META.PublicRelease.300517.nodup.QC.txt --target PAPG21-3_capstone_1000g_allphase3SAS_GRCh37.QC-merged_excludemissnp.final --no-default --chr chr --A1 effect_allele --A2 noneffect_allele --snp snptestid --pvalue p-value_gc --bp bp_hg19 --or --cov PAPG21-3_capstone1000G_merged.QC.COV.final  --pheno PAPG21-3_capstone_1000g_allphase3SAS_GRCh37.QC-merged_excludemissnp.pheno.final --pheno-col T2D --ignore-fid --missing CENTER --thread 2  --stat logOR   --binary-target T   --score avg   --quantile 100  --quant-break 1,5,10,20,40,60,80,90,95,99,100   --quant-ref 60  --out UKBB_EXOME_CAD_PRSice_avg_PAPG21-3_1000gSAS_GRCh37
