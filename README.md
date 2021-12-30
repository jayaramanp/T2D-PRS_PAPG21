
# 1. Base Data


####T2D
#### download BASE GWAS summary stats from DIAGRAM (http://diagram-consortium.org/downloads.html) for T2D from Mahajan et al 2018
unzip T2D.TranEthnic.zip

## used the BMI adjusted.txt

###### DO NOT LIFTOVER TO 38, USE ORIGINAL GrCh37 since we have lifted over our target data to 37. 
#### it was found to be in GrCh37 assembly. hence we needed to perform liftover and also ensure the GWAS summary stats format remained the same. 

#### used this: https://github.com/hakyimlab/summary-gwas-imputation/wiki/GWAS-Harmonization-And-Imputation

#module load python/3.8.2

#python summary-gwas-imputation/src/gwas_parsing.py -gwas_file T2D_TranEthnic.BMIadjusted.txt -output_column_map rsID variant_id -output_column_map OTHER_ALLELE non_effect_allele -output_column_map EFFECT_ALLELE effect_allele -output_column_map BETA effect_size -output_column_map Pvalue pvalue -output_column_map SE standard_error -output_column_map SNP chr_pos --chromosome_format -split_column chr_pos ':' chromosome position -output_column_map Neff sample_size -output_order variant_id panel_variant_id chromosome position effect_allele non_effect_allele frequency pvalue zscore effect_size standard_error sample_size n_cases -liftover hg19ToHg38.over.chain.gz -output T2D_TranEthnic.BMIadjusted.liftoverhg38.txt
### GrCh38 Base data: T2D_TranEthnic.BMIadjusted.liftoverhg38.txt



###### QC on GrCh37 for the Base data

## remove dups from Base data
cat T2D_TranEthnic.BMIadjusted.txt | awk '{seen[$2]++; if(seen[$2]==1){ print}}' >> T2D_TranEthnic.BMIadjusted.nodup.txt


### remove ambiguous snps:
cat T2D_TranEthnic.BMIadjusted.nodup.txt | awk '!( ($5=="A" && $6=="T") || ($5=="T" && $6=="A") || ($5=="G" && $6=="C") || ($5=="C" && $6=="G")) {print}' > T2D_TranEthnic.BMIadjusted.nodup.QC.txt



### remove first column from file - its got no data anyway!
awk '{$1=""}1' T2D_TranEthnic.BMIadjusted.nodup.QC.txt >> T2D_TranEthnic.BMIadjusted.nodup.QC.final.txt


#now its ready to be used!


### B. MAGIC dataset

##removedup snps
cat MAGIC1000G_FG_SAS.tsv | awk '{seen[$1]++; if(seen[$1]==1){ print}}' >> MAGIC1000G_FG_SAS.nodup.txt

### remove ambiguous snps
cat MAGIC1000G_FG_SAS.nodup.txt | awk '!( ($4=="A" && $5=="T") || ($4=="T" && $5=="A") || ($4=="G" && $5=="C") || ($4=="C" && $5=="G")) {print}' > MAGIC1000G_FG_SAS.nodup.QC.txt

### instead of renaming headers. maybe add additional params in your PRSice command



#### C. METAANALYSIS DIAGRAM SCott et al 2017

## split first column into chr and bp and keep snp column as well. 

awk '{split($1,s,":"); print s[1],s[2],$1,$2,$3,$4,$5,$6,$7;}' METAANALYSIS_DIAGRAM_SE1.BMI.txt > METAANALYSIS_DIAGRAM_SE1.BMI.MOD.txt

## make file tab delimited
use vi and %s/\s/\t/g

## column has effect - which is same as beta

### remove duplicated snps
cat METAANALYSIS_DIAGRAM_SE1.BMI.MOD.txt | awk '{seen[$3]++; if(seen[$3]==1){ print}}' >> METAANALYSIS_DIAGRAM_SE1.BMI.MOD.nodup.txt


### remove ambigous snps
cat METAANALYSIS_DIAGRAM_SE1.BMI.MOD.nodup.txt | awk '!( ($4=="A" && $5=="T") || ($4=="T" && $5=="A") || ($4=="G" && $5=="C") || ($4=="C" && $5=="G")) {print}' > METAANALYSIS_DIAGRAM_SE1.BMI.MOD.nodup.QC.txt


######## D. UKBB meta analysis A meta-analysis of UK Biobank SOFT CAD GWAS (interim release) with CARDIoGRAMplusC4D 1000 Genomes-based GWAS and the Myocardial Infarction Genetics and CARDIoGRAM Exome

###remove dups:
cat UKBB.GWAS1KG.EXOME.CAD.SOFT.META.PublicRelease.300517.txt | awk '{seen[$2]++; if(seen[$2]==1){ print}}' >> UKBB.GWAS1KG.EXOME.CAD.SOFT.META.PublicRelease.300517.nodup.txt

### remove ambigous snps:
cat UKBB.GWAS1KG.EXOME.CAD.SOFT.META.PublicRelease.300517.nodup.txt | awk '!( ($5=="A" && $6=="T") || ($5=="T" && $6=="A") || ($5=="G" && $6=="C") || ($5=="C" && $6=="G")) {print}' > UKBB.GWAS1KG.EXOME.CAD.SOFT.META.PublicRelease.300517.nodup.QC.txt


##### manually modify snp id position for
Error: Invalid loci for rs188690587: 9.3e+07


#### turns out there are many more:
grep 'e+0' UKBB.GWAS1KG.EXOME.CAD.SOFT.META.PublicRelease.300517.nodup.QC.txt 
2:72000000_G_A	rs61573637	2	7.2e+07	G	A	0.88523	0.0065	0.013	0.616807	330346	no	0.98
6:136000000_T_C	rs55926138	6	1.36e+08	C	T	0.0188	0.01293	0.03557	0.716188	293919	no	0.98
12:48000000_T_C	rs855274	12	4.8e+07	C	T	0.88228	0.01649	0.01407	0.2411	334732	no	1
14:81000000_TC_T	rs151019276	14	8.1e+07	TC	T	0.96264	0.03614	0.03261	0.267785	255037	no	0.66


###### E. UACR DM cohort:

###remove rows that dont have RS ID
cat UACR_formatted_20170020-UACR_DM-All-nstud_18-SumMac_400.tbl.rsid | awk 'NF==10{print}{}' >> UACR_formatted_20170020_DM-All-nstud_18-SumMac_400.tbl.ALLrsid


### remove dups
cat UACR_formatted_20170020_DM-All-nstud_18-SumMac_400.tbl.ALLrsid | awk '{seen[$3]++; if(seen[$3]==1){ print}}' >> UACR_formatted_20170020_DM-All-nstud_18-SumMac_400.tbl.ALLrsid.nodup.txt

##remove ambiguous snps
cat UACR_formatted_20170020_DM-All-nstud_18-SumMac_400.tbl.ALLrsid.nodup.txt | awk '!( ($4=="A" && $5=="T") || ($4=="T" && $5=="A") || ($4=="G" && $5=="C") || ($4=="C" && $5=="G")) {print}' > UACR_formatted_20170020_DM-All-nstud_18-SumMac_400.tbl.ALLrsid.nodup.QC.txt



### 2. Target data processing starts here:

#### FINAL update dont use Crossmap!!
#also liftover VCF file to GrCh37 hg19 using CRossMap
### http://crossmap.sourceforge.net/#chain-file

module load python/3.7.3
module load crossmap/0.2.9

gzip -d PAPG21-3.vcf.gz

CrossMap.py vcf hg38ToHg19.over.chain.gz PAPG21-3.vcf ../reference/Homo_sapiens_assembly38.fasta PAPG21-3_capstone.vcf


###FINAL update - this works!
###### also tried PICARD liftOver:
module load picard

java -jar $PICARD LiftoverVcf I=PAPG21-3.vcf O=PAPG21-3_liftover2hg19.vcf CHAIN=hg38ToHg19.over.chain.gz REJECT=PAPG21-3_liftover2hg19_rejected_variants.vcf R=../reference/GRCh37.p13.genome.fa WARN_ON_MISSING_CONTIG=true


##### sort VCF - if PICARD  no need to sort VCF

grep "^#" PAPG21-3_capstone.vcf > PAPG21-3_capstone.sorted.vcf && grep -v "^#" PAPG21-3_capstone.vcf | sort -V -k1,1 -k2,2n >> PAPG21-3_capstone.sorted.vcf



### THIS WORKS!

##### create PLINK binary files:

plink2 --vcf PAPG21-3_liftover2hg19.vcf.gz --max-alleles 2 --rm-dup --make-bed --out PAPG21-3_capstone

you will see it doesnt find sex = female. 


####change fam file to this:
echo "0"$'\t'"PAPG21-3_capstone"$'\t'"0"$'\t'"0"$'\t'"2"$'\t'"2"$'\n' > PAPG21-3_capstone.fam

###rerun plink2
plink2 --bfile PAPG21-3_capstone --max-alleles 2 --rm-dup --make-bed --out PAPG21-3_capstone

## remove dup snps:

cut -f 2 PAPG21-3_capstone.bim | sort | uniq -d > snp.dups

## rerun plink2
plink2 --bfile PAPG21-3_capstone --max-alleles 2 --allow-extra-chr --exclude snp.dups --rm-dup --make-bed --out PAPG21-3_capstone_nodup


### QC using plink for patient sample:
plink --bfile PAPG21-3_capstone_nodup --maf 0.01 --mind 0.1 --geno 0.1 --make-bed --out PAPG21-3_capstone_nodup.qc


#### this is now ready to use. you will need 1000G data to add power to this. 

# 2. Target Data

ml plink2
ml plink


#### 1000G data

## you can download 1000G vcf and then use the commands here to create plink files: https://www.cog-genomics.org/plink2/data#merge3
## or easier:
# download 1000G Phase3 binary files data for all super populations from Plink2 website: from here: https://www.cog-genomics.org/plink/2.0/resources#1kg_phase3
## download the bold links. Unzip the file like the command above in the url
 
#pull out list of samples you need to keep then run: 
plink2 --pfile all_phase3 --max-alleles 2 --make-bed --keep keep_SAS.list --out all_phase3_SAS


##NOW run QC. #Standard GWAS QC  
plink --bfile all_phase3_SAS --allow-extra-chr --maf 0.01 --hwe 1e-6 --geno 0.01 --mind 0.01 --write-snplist --make-just-fam --out all_phase3_SAS.QC

## prune highly correlated SNPs
plink --bfile all_phase3_SAS --allow-extra-chr --keep all_phase3_SAS.QC.fam --extract all_phase3_SAS.QC.snplist --indep-pairwise 200 50 0.25 --out all_phase3_SAS.QC

## compute heterozygosity rates
plink --bfile all_phase3_SAS --allow-extra-chr --extract all_phase3_SAS.QC.prune.in --keep all_phase3_SAS.QC.fam --het --out all_phase3_SAS.QC



### the files you see so far:
-rw-r----- 1 jayarp02 bsr2401  54K Dec 24 22:29 all_phase3.psam
-rw-r----- 1 jayarp02 bsr2401 2.3G Dec 24 22:29 all_phase3.pgen.zst?dl=1
-rw-r----- 1 jayarp02 bsr2401 1.3G Dec 24 22:30 all_phase3.pvar.zst?dl=1
-rw-r----- 1 jayarp02 bsr2401 6.3G Dec 24 22:31 all_phase3.pgen
-rw-r----- 1 jayarp02 bsr2401  12G Dec 24 22:32 all_phase3.pvar
-rw-r----- 1 jayarp02 bsr2401  11K Dec 24 22:36 all_SAS_1KG_phase3.psam
-rw-r----- 1 jayarp02 bsr2401 3.9K Dec 24 22:53 keep_SAS.list
-rw-r----- 1 jayarp02 bsr2401 9.1K Dec 24 22:54 all_phase3_SAS.fam
-rw-r----- 1 jayarp02 bsr2401 2.4G Dec 24 22:54 all_phase3_SAS.bim
-rw-r----- 1 jayarp02 bsr2401 9.6G Dec 24 22:54 all_phase3_SAS.bed
-rw-r----- 1 jayarp02 bsr2401 1023 Dec 24 22:54 all_phase3_SAS.log
-rw-r----- 1 jayarp02 bsr2401 106M Dec 24 23:13 all_phase3_SAS.QC.snplist
-rw-r----- 1 jayarp02 bsr2401 9.1K Dec 24 23:13 all_phase3_SAS.QC.fam
-rw-r----- 1 jayarp02 bsr2401  14M Dec 24 23:20 all_phase3_SAS.QC.prune.in
-rw-r----- 1 jayarp02 bsr2401  96M Dec 24 23:20 all_phase3_SAS.QC.prune.out
drwxr-s--- 2 jayarp02 bsr2401 4.0K Dec 24 23:22 .
-rw-r----- 1 jayarp02 bsr2401  32K Dec 24 23:22 all_phase3_SAS.QC.het
-rw-r----- 1 jayarp02 bsr2401 1.2K Dec 24 23:22 all_phase3_SAS.QC.log

### remove F coeff more than 3SD from mean
Rscript remove_fcoeff_3SD_frommean.R 


### sexcheck
plink --bfile all_phase3_SAS --allow-extra-chr --extract all_phase3_SAS.QC.prune.in --keep all_phase3_SAS.valid.sample --check-sex --out all_phase3_SAS.QC

### if the file .sexcheck has PROBLEM then remove that from file by reruning plink with remove

plink --bfile all_phase3_SAS --allow-extra-chr --extract all_phase3_SAS.QC.prune.in --keep all_phase3_SAS.valid.sample_removedSexcheckQC --check-sex --out all_phase3_SAS.QC

### create final valid file:

Rscript create_FID_IID_validsamples.R

## then rerun for Relatedness:
plink --bfile all_phase3_SAS --allow-extra-chr --extract all_phase3_SAS.QC.prune.in --keep all_phase3_SAS.valid.sample_removedSexcheckQC --rel-cutoff 0.125 --out all_phase3_SAS.QC



###SIMULATE PHENOTYPE
## The files above have no phenotype. simulate 1000G phenotypes via this link: https://github.com/bambrozio/bioinformatics/blob/master/utils/1k-genomes-phenotype-simulator.ipynb

### download panel file:
wget ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/integrated_call_samples_v3.20130502.ALL.panel

###extract SAS super pop sample rows. 

###R script to simulate phenotypes for SAS:
ml R
R create_phenotypes.R 

### then merge the simulated phenotype (1=ctrl, 2=case) into the fam file

Rscript pull_simPhenotype_into_fam.R

####This will also create a new covariates file and a Phenotype file for use later!

####copy the all_phase3_SAS.QC.fam to all_phase3_SAS.fam (after backing up the original file with _orig)

####then generate final set of QCs with Phenotype:

plink --bfile all_phase3_SAS --allow-extra-chr --make-bed --keep all_phase3_SAS.QC.rel.id --out all_phase3_SAS.QC.final --extract all_phase3_SAS.QC.snplist


####MERGE:
plink --bfile 1000G/all_phase3_SAS.QC.final --allow-extra-chr --bmerge PAPG21-3_capstone_nodup.qc --make-bed --out PAPG21-3_capstone_1000g_allphase3SAS_GRCh37.QC

###you may see merge errors with variants seen before. refer to this: https://www.cog-genomics.org/plink2/data#merge3
##step1:
plink --bfile PAPG21-3_capstone_nodup.qc --flip PAPG21-3_capstone_1000g_allphase3SAS_GRCh37.QC-merge.missnp --make-bed --out PAPG21-3_capstone_nodup.qc_trial

##step2:
plink --bfile 1000G/all_phase3_SAS.QC.final --allow-extra-chr --bmerge PAPG21-3_capstone_nodup.qc_trial --make-bed --out PAPG21-3_capstone_1000g_allphase3SAS_GRCh37_merged_trial

### THIS did not resolve the issue!



#### alt approach because the above approach doesnt work! exclude the snps
plink --bfile 1000G/all_phase3_SAS.QC.final --allow-extra-chr --exclude merge_test/PAPG21-3_capstone_1000g_allphase3SAS_GRCh37.QC-merge.missnp --make-bed --out merge_test/all_phase3_SAS.QC.final_tmp

plink --bfile PAPG21-3_capstone_nodup.qc --allow-extra-chr --exclude merge_test/PAPG21-3_capstone_1000g_allphase3SAS_GRCh37.QC-merge.missnp --make-bed --out merge_test/PAPG21-3_capstone_nodup.qc_tmp

plink --bfile merge_test/all_phase3_SAS.QC.final_tmp --allow-extra-chr --bmerge merge_test/PAPG21-3_capstone_nodup.qc_tmp --make-bed --out PAPG21-3_capstone_1000g_allphase3SAS_GRCh37.QC-merged_excludemissnp

## you can then remove the _tmp files
rm merge_test/all_phase3_SAS.QC.final_tmp.*
rm merge_test/PAPG21-3_capstone_nodup.qc_tmp.*


#METHOD1 generate new PRS:


###calculate eigenvalues
plink --bfile 1000G/all_phase3_SAS.QC.final --allow-extra-chr --pca --out 1000G/all_phase3_SAS.QC.final.pca


#### create COVARIATE file - not right! the R script above creates the covariates and Phenotypes file. 
#cat 1000G/all_phase3_SAS.QC.final.fam | cut -d ' ' -f 1-2,6 > 1000G/all_phase3_SAS.QC.final.covariates

## R file to plot PCA eigenvectors for SAS file only
Rscript plot_PCA_allphase3SASQCFinal.R 


#### calculate eigen values for PAPG_1000G merged 
plink --bfile PAPG21-3_capstone_1000g_allphase3SAS_GRCh37.QC-merged_excludemissnp --allow-extra-chr --pca --out PAPG21-3_capstone_1000g_allphase3SAS_GRCh37.QC-merged_excludemissnp.pca

##### create COVARIATE file for PAPG_1000G_merged - again not right! will rename this to Phenotype and manually add random values for covariates. Copy covariates from all_phase3_SAS.QC.NEWCOVARIATES and manually add one extra row for the patient sample.  
cat PAPG21-3_capstone_1000g_allphase3SAS_GRCh37.QC-merged_excludemissnp.fam | cut -d ' ' -f 1-2,6 > PAPG21-3_capstone_1000g_allphase3SAS_GRCh37.QC-merged_excludemissnp.Pheno

#### R file to plot PCA eigenvectors for PAPG_1000G_merged 
Rscript plot_PCA_PAPG_merged_1000GSASFinal.R 


##### run PRSice

## but first remove any unnatural extra chr:
plink --bfile PAPG21-3_capstone_1000g_allphase3SAS_GRCh37.QC-merged_excludemissnp --allow-extra-chr --not-chr PAR1 --make-bed --out PAPG21-3_capstone_1000g_allphase3SAS_GRCh37.QC-merged_excludemissnp.final

### also remove first column in pheno file because no FID is present. 
awk '{$1=""}1' PAPG21-3_capstone_1000g_allphase3SAS_GRCh37.QC-merged_excludemissnp.pheno >> PAPG21-3_capstone_1000g_allphase3SAS_GRCh37.QC-merged_excludemissnp.pheno.final

#### maybe also remove from COV file
awk '{$1=""}1' PAPG21-3_capstone1000G_merged.QC.COV >> PAPG21-3_capstone1000G_merged.QC.COV.final



##### 1. with DIAGRAM 2018 T2D Mahajan et al. 
runPRSice_T2D.sh

#### MiSIGDB for pathways and gene sets
 bsub< LSF_plotT2D_MiSIGDB_submit

### uses the script: plot_T2D_MiSIGDB.sh



#### 2. with MAGIC
runPRSice_MAGIC.sh
Rscript plotPRS_MAGIC.R

##### MiSIGDB for pathways and gene sets:
sh plot_MAGIC_MiSIGDB.sh

### but the above is slow (LSF is a better strategy)




#### 3. with DIAGRAM 2017 data Scott et al.
runPRSice_DIAGRAM2017.sh
Rscript plotPRS_DIAGRAM2017.R





##### 4. with UKBB EXome CAD
runPRSice_UKBB_CAD.sh
Rscript plotPRS_UKBB.R

##### 5. run with UACR DM
runPRSice_UACRDM.sh
Rscript plotPRS_UACRDM.R



##### plot all PRS:
plot_all_MAGICT2DDIAGRAMUKBBUACR_PRS.R




#Method2
#### PRS catalogues file: https://www.pgscatalog.org/score/PGS000036/

compare against PGS000036 and PGS000804

