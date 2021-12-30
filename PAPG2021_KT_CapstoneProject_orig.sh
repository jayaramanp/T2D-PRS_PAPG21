# PAPG2021 Capstone Project 
# auth: Kayla Townsley 

# PRSice to calculate trait polygenic risk score and PRSet to caluclate Atopic Derm/BIP pathway-based polygenic risk for immunological signature gene sets (C2 from MSigDB)

#1000genomes EUR reference genome
scp -r /Users/townsk05/Desktop/EUR /sc/arion/projects/bsr*/PA*/townsk05_dir/PRsice/source_data

cd /sc/arion/projects/bsr*/PA*/townsk05_dir/PRsice/source_data

#QC of Base Data (GWAS Summary Stats)

##Heritability check : h2snp > 0.05 by LDSC
##Effect Allele: importat to know which allele is effect allele for PRS association results to be in the correct direction
##Genome build
##Standard GWAS QC

gunzip -c /sc/arion/projects/psychgen/pgc/psychchip/finalRelease/newpreimputation/bip57/distribution/pgc3_bip/daner_bip_pgc3.gz|\
awk 'NR==1 || ($6 > 0.01) && ($8 > 0.8) {print}' |\
gzip  > daner_bip_pgc3.gz

##Remove duplicate SNPS

gunzip -c daner_bip_pgc3.gz |\
awk '{print $2}' |\
sort |\
uniq -d > duplicated.snp

#no duplicated SNP IDs
#gunzip -c daner_bip_pgc3.gz |\
#grep -vf duplicated.snp |\
#gzip - > daner_bip_pgc3.nodup.gz

#Remove ambiguous SNPS 

gunzip -c daner_bip_pgc3.gz |\
awk '!( ($4=="A" && $5=="T") || \
        ($4=="T" && $5=="A") || \
        ($4=="G" && $5=="C") || \
        ($4=="C" && $5=="G")) {print}' |\
    gzip > daner_bip_pgc3.QC.gz

#1087297 ambiguos snps removed (6085542 SNPs remaining in QC'd file)

#Liftover of personal genome VCF from Builg hg38 to build hg19
CrossMap.py vcf GRCh38_to_GRCh37.chain PAPG-townsk05*.vcf hg19.fa PAPG-townsk05_GRCh37

#create bed files from personal genome vcfs
ml plink 
ml plink2
plink2 --vcf PAPG-townsk05_GRCh37 --max-alleles 2 --rm-dup --make-bed --out PAPG-townsk05_GRCh37
cut -f 2 PAPG-townsk05_GRCh37.bim | sort | uniq -d > snp.dups
plink --bfile PAPG-townsk05_GRCh37 --allow-extra-chr --exclude snp.dups --make-bed --out PAPG-townsk05_GRCh37_nodups


# Quality Control of Personal Genome

plink --bfile PAPG-townsk05_GRCh37_nodups \
    --maf 0.01 \
    --mind 0.1 \
    --geno 0.1 \
    --make-bed \
    --out PAPG-townsk05_GRCh37.qc

#QC Target Data
    # Sample size should be > 100
    # Genome build must match between target and base data
#Standard GWAS QC 
ml plink
plink \
    --bfile EUR \
    --maf 0.01 \
    --hwe 1e-6 \
    --geno 0.01 \
    --mind 0.01 \
    --write-snplist \
    --make-just-fam \
    --out EUR.QC

 # prune highly correlated SNPs
 plink \
    --bfile EUR \
    --keep EUR.QC.fam \
    --extract EUR.QC.snplist \
    --indep-pairwise 200 50 0.25 \
    --out EUR.QC

 # compute heterozygosity rates
 plink \
    --bfile EUR \
    --extract EUR.QC.prune.in \
    --keep EUR.QC.fam \
    --het \
    --out EUR.QC

 # Remove individuals with F coefficients that are more than 3 SDs from the mean
 ml R
 R
library(data.table)
# Read in file
dat <- fread("EUR.QC.het")
# Get samples with F coefficient within 3 SD of the population mean
valid <- dat[F<=mean(F)+3*sd(F) & F>=mean(F)-3*sd(F)] 
# print FID and IID for valid samples
fwrite(valid[,c("FID","IID")], "EUR.valid.sample", sep="\t") 
q() # exit R

# Mismatching SNPs
## magrittr allow us to do piping, which help to reduce the 
# amount of intermediate data types
#library(data.table)
#library(magrittr)
# Read in bim file 
#bim <- fread("PAPG-townsk05_1000gEUR_hg19.qc2.bim") %>%
    # Note: . represents the output from previous step
    # The syntax here means, setnames of the data read from
    # the bim file, and replace the original column names by 
    # the new names
    #setnames(., colnames(.), c("CHR", "SNP", "CM", "BP", "B.A1", "B.A2")) %>%
    # And immediately change the alleles to upper cases
    #.[,c("B.A1","B.A2"):=list(toupper(B.A1), toupper(B.A2))]
# Read in summary statistic data (require data.table v1.12.0+)
 # And immediately change the alleles to upper cases
 # Read in QCed SNPs
#BIP <- fread("../data/daner_bip_pgc3.QC.gz") %>% .[,c("A1","A2"):=list(toupper(A1), toupper(A2))]
#qc <- fread("EUR.QC.snplist", header=F)

# Merge summary statistic with target
# And filter out QCed SNPs
#info <- merge(bim, BIP, by=c("SNP", "CHR", "BP")) %>%
 #   .[SNP %in% qc[,V1]]

# Function for calculating the complementary allele
#complement <- function(x){
 #   switch (x,
  #      "A" = "T",
   #     "C" = "G",
    #    "T" = "A",
     #   "G" = "C",
      #  return(NA)
#    )
#} 
# Get SNPs that have the same alleles across base and target
#info.match <- info[A1 == B.A1 & A2 == B.A2, SNP]
# Identify SNPs that are complementary between base and target
#com.snps <- info[sapply(B.A1, complement) == A1 &
                    sapply(B.A2, complement) == A2, SNP]
# Now update the bim file
#bim[SNP %in% com.snps, c("B.A1", "B.A2") :=
 #       list(sapply(B.A1, complement),
  #          sapply(B.A2, complement))]
# identify SNPs that need recoding
#recode.snps <- info[B.A1==A2 & B.A2==A1, SNP]
# Update the bim file
#bim[SNP %in% recode.snps, c("B.A1", "B.A2") :=
        list(B.A2, B.A1)]

# identify SNPs that need recoding & complement
#com.recode <- info[sapply(B.A1, complement) == A2 &
                    sapply(B.A2, complement) == A1, SNP]
# Now update the bim file
#bim[SNP %in% com.recode, c("B.A1", "B.A2") :=
        list(sapply(B.A2, complement),
            sapply(B.A1, complement))]
# Write the updated bim file
fwrite(bim[,c("SNP", "B.A1")], "EUR.a1", col.names=F, sep="\t")

mismatch <- bim[!(SNP %in% info.match |
                    SNP %in% com.snps |
                    SNP %in% recode.snps |
                    SNP %in% com.recode), SNP]
write.table(mismatch, "EUR.mismatch", quote=F, row.names=F, col.names=F)
q() # exit R

# Sex chromosomes 
plink \
    --bfile EUR \
    --extract EUR.QC.prune.in \
    --keep EUR.valid.sample \
    --check-sex \
    --out EUR.QC
R
# Read in file
valid <- read.table("EUR.valid.sample", header=T)
dat <- read.table("EUR.QC.sexcheck", header=T)
valid <- subset(dat, STATUS=="OK" & FID %in% valid$FID)
write.table(valid[,c("FID", "IID")], "EUR.QC.valid", row.names=F, col.names=F, sep="\t", quote=F) 
q() # exit R

#Relatedness

plink \
    --bfile EUR \
    --extract EUR.QC.prune.in \
    --keep EUR.QC.valid \
    --rel-cutoff 0.125 \
    --out EUR.QC

#Generate final QC'ed target data file
plink \
    --bfile EUR \
    --make-bed \
    --keep EUR.QC.rel.id \
    --out EUR.QC \
    --extract EUR.QC.snplist \

# Merge QC' EUR reference genotypes with personal genome

plink --bfile ../EUR/EUR.QC --bmerge PAPG-townsk05_GRCh37.qc --make-bed --out PAPG-townsk05_1000gEUR_GRCh37.QC

#You may need to exclude multiallelic variants
plink --bfile PAPG-townsk05_GRCh37.qc --exclude  PAPG-townsk05_1000gEUR_GRCh37.QC-merge.missnp --make-bed --out PAPG-townsk05_GRCh37.qc


#Calculate eigenvalues (PCA) 
plink --bfile EUR.QC --pca --out EUR.QC.pca

# Plot PCAs
ml R
R 
eur.pca<- read.table("EUR.QC.pca.eigenvec")
eur.pca <- as.data.frame(eur.pca[3:22])
colnames(eur.pca) <- c("PCA_1", "PCA_2", "PCA_3", "PCA_4", "PCA_5", "PCA_6", "PCA_7", "PCA_8", "PCA_9", "PCA_10", "PCA_11", "PCA_12", "PCA_13", "PCA_14", "PCA_15", "PCA_16", "PCA_17", "PCA_18", "PCA_19", "PCA_20")

library(ggfortify) 
pca_res <- prcomp(eur.pca,  scale = TRUE) 
summary(pca_res)
var_explained <-pca_res$sdev^2/sum(pca_res$sdev^2)
var_explained[1:5]
PCA <- as.data.frame(pca_res$x)

tiff("EUR_QC_pca.tiff", units="in", width=10, height=10, res=300)
ggplot(PCA, aes(x=PC1,y=PC2)) + geom_point(size=4) + theme_bw(base_size=24) + labs(x=paste0("PC1: ",round(var_explained[3]*100,1),"%"), y=paste0("PC2: ",round(var_explained[4]*100,1),"%")) + theme(legend.position="top")
dev.off()

covariate <- read.table("EUR.cov", header=T)
pcs <- read.table("EUR.QC.pca.eigenvec", header=F)
colnames(pcs) <- c("FID","IID", paste0("PC",1:6))
cov <- merge(covariate, pcs, by=c("FID", "IID"))
write.table(cov,"EUR.covariate", quote=F, row.names=F)
q()

#Calculate eigenvalues (PCA) 
plink --bfile PAPG-townsk05_1000gEUR_GRCh37.QC --pca --out PAPG-townsk05_1000gEUR_GRCh37.QC.pca

# Plot PCAs
ml R
R 
pca_TargetData <- read.table("PAPG-townsk05_1000gEUR_hg19.QC.pca.eigenvec")
TD_pca <- as.data.frame(pca_TargetData[3:22])
colnames(TD_pca) <- c("PCA_1", "PCA_2", "PCA_3", "PCA_4", "PCA_5", "PCA_6", "PCA_7", "PCA_8", "PCA_9", "PCA_10", "PCA_11", "PCA_12", "PCA_13", "PCA_14", "PCA_15", "PCA_16", "PCA_17", "PCA_18", "PCA_19", "PCA_20")

library(ggfortify) 
pca_res <- prcomp(TD_pca,  scale = TRUE) 
summary(pca_res)
var_explained <-pca_res$sdev^2/sum(pca_res$sdev^2)
var_explained[1:5]
PCA <- as.data.frame(pca_res$x)

tiff("TargetData_pca.tiff", units="in", width=10, height=10, res=300)
ggplot(PCA, aes(x=PC1,y=PC2)) + geom_point(color = "blue", size=4) + theme_bw(base_size=24) + labs(x=paste0("PC1: ",round(var_explained[3]*100,1),"%"), y=paste0("PC2: ",round(var_explained[4]*100,1),"%")) + theme(legend.position="top")
dev.off()

covariate <- read.table("PAPG-05_100gEUR.cov", header=T)
pcs <- read.table("PAPG-townsk05_1000gEUR_GRCh37.QC.pca.eigenvec", header=F)
colnames(pcs) <- c("FID","IID", paste0("PC",1:20))
cov <- merge(covariate, pcs, by=c("FID", "IID"))
write.table(cov,"PAPG-townsk05_1000gEUR_GRCh37.QC.covariate", quote=F, row.names=F)
q()

# Calcualting PRS score for Binary traits from Personal + Reference Genomes
Rscript ../../PRSice.R --dir .\
    --prsice ../../PRSice_linux \
    --base  ../../data/AtopicDerm.sumstats \
    --target ../../data/PAPG-townsk05_1000gEUR_GRCh37.QC \
    --cov ../../data/PAPG-05_1000gEUR.cov  \
    --pheno ../../data/pheno \
    --ignore-fid \
    --thread 1 \
    --stat BETA \
    --binary-target T \
    --score avg \
    --quantile 100 \
    --quant-break 1,5,10,20,40,60,80,90,95,99,100 \
    --quant-ref 60 \
    --out AtopicDerm_PRSice_avg_PAPG-townsk05_1000gEUR_GRCh37

    Rscript ../../PRSice_linux/PRSice.R --dir .\
    --prsice ../../PRSice_linux/PRSice_linux \
    --base  ../../data/AtopicDerm.sumstats \
    --target ../../EUR/EUR.QC \
    --pheno ../../EUR/EUR.QC.pheno \
    --ignore-fid \
    --thread 1 \
    --stat BETA \
    --binary-target T \
    --score avg \
    --quantile 100 \
    --quant-break 1,5,10,20,40,60,80,90,95,99,100 \
    --quant-ref 60 \
    --out AtopicDerm_PRSice_avg_1000gEUR_GRCh37


    Rscript ../../PRSice_linux/PRSice.R --dir .\
    --prsice ../../PRSice_linux/PRSice_linux \
    --base  ../../data/daner_bip_pgc3.QC.txt \
    --target ../../EUR/EUR.QC \
    --pheno ../../EUR/EUR.QC.pheno \
    --ignore-fid \
    --thread 1 \
    --stat OR \
    --binary-target T \
    --score con-std \
    --quantile 100 \
    --quant-break 1,5,10,20,40,60,80,90,95,99,100 \
    --quant-ref 60 \
    --out PRSice_con-std_1000gEUR_hg19_BIP3

# Create a density plot of the PRS score distribution

ml R
R
prs.all <- read.table("PRSice_cont-std_PAPG-townsk05_1000gEUR_hg19_BIP3.best", header = TRUE)
prs.all <- as.data.frame(prs.all)
# Find the mean
library(plyr)
library(ggplot2)
cprs <- (PRS.mean=mean(prs.all$PRS))
cprs

tiff("BIP3_PRS_densityplot.tiff", units="in", width=10, height=10, res=300)
ggplot(prs.all, aes(x=PRS)) + geom_density() + geom_vline(data=prs.all, aes(xintercept=PRS.mean), linetype="dashed", size=1) + geom_vline(data=prs.all, aes(xintercept=-0.418322903), color = "blue", size = 0.5)
dev.off()



#PRSet
# Msigdb C2: CGP (Chemical and genetic perturbations; n=3358 gene sets)
 Rscript ../../PRSice_linux/PRSice.R --dir .\
    --prsice ../../PRSice_linux/PRSice_linux \
    --base ../../data/daner_bip_pgc3.QC.txt  \
    --target ../../data/PAPG-townsk05_1000gEUR_hg19.QC \
    --pheno ../../data/PAPG-townsk05_1000gEUR.BIP.pheno \
    --ignore-fid \
    --binary-target T \
    --thread 1 \
    --gtf ../../source_data/Homo_sapiens.GRCh37.87.gtf \
    --msigdb ../../source_data/c2.cgp.v7.2.symbols.gmt \
    --multi-plot 10 \
    --out PRSet_BIP3_PAPG-townsk05_1000g_EUR_hg19_c2.cgp


# Msigdb C2: All curated gene sets (CGP & Canonical Pathways (CP) from BIOCARTA, KEGG, PID, REACTOME,and WikiPathways; n=6226 gene sets)
 Rscript ../../PRSice_linux/PRSice.R --dir .\
    --prsice ../../PRSice_linux/PRSice_linux \
    --base ../../data/daner_bip_pgc3.QC.txt  \
    --target ../../EUR/EUR.QC \
    --pheno ../../EUR/EUR.QC.pheno \
    --ignore-fid \
    --binary-target T \
    --thread 1 \
    --gtf ../../source_data/Homo_sapiens.GRCh37.87.gtf \
    --msigdb ../../source_data/c2.all.v7.2.symbols.gmt \
    --multi-plot 10 \
    --out PRSet_BIP3_EUR_hg19_c2.all


# Msigdb C3: Regulatory target gene sets; n=3556 gene sets)
 Rscript ../../PRSice_linux/PRSice.R --dir .\
    --prsice ../../PRSice_linux/PRSice_linux \
    --base ../../data/daner_bip_pgc3.QC.txt  \
    --target ../../EUR/EUR.QC \
    --pheno ../../EUR/EUR.QC.pheno \
    --ignore-fid \
    --binary-target T \
    --thread 1 \
    --gtf ../../source_data/Homo_sapiens.GRCh37.87.gtf \
    --msigdb ../../source_data/c3.all.v7.2.symbols.gmt \
    --multi-plot 10 \
    --out PRSet_BIP3_EUR_hg19_c3.regulatory-targets

# Msigdb C5: Ontology gene sets ; n=14765 gene sets)
 Rscript ../../PRSice_linux/PRSice.R --dir .\
    --prsice ../../PRSice_linux/PRSice_linux \
    --base ../../data/daner_bip_pgc3.QC.txt  \
    --target ../../EUR/EUR.QC \
    --pheno ../../EUR/EUR.QC.pheno \
    --ignore-fid \
    --binary-target T \
    --thread 1 \
    --gtf ../../source_data/Homo_sapiens.GRCh37.87.gtf \
    --msigdb ../../source_data/c5.all.v7.2.symbols.gmt \
    --multi-plot 10 \
    --out PRSet_BIP3_EUR_hg19_c5.GO

# Msigdb C7: Immunologic signature gene sets; n=4872 gene sets)
 Rscript ../../PRSice_linux/PRSice.R --dir .\
    --prsice ../../PRSice_linux/PRSice_linux \
    --base ../../data/daner_bip_pgc3.QC.txt  \
    --target ../../EUR/EUR.QC \
    --pheno ../../EUR/EUR.QC.pheno \
    --ignore-fid \
    --binary-target T \
    --thread 1 \
    --gtf ../../source_data/Homo_sapiens.GRCh37.87.gtf \
    --msigdb ../../source_data/c7.all.v7.2.symbols.gmt \
    --multi-plot 10 \
    --out PRSet_BIP3_EUR_hg19_c7.immune

# Msigdb H: Hallmark gene sets; n=50 gene sets)
 Rscript ../../PRSice_linux/PRSice.R --dir .\
    --prsice ../../PRSice_linux/PRSice_linux \
    --base ../../data/daner_bip_pgc3.QC.txt  \
    --target ../../EUR/EUR.QC \
    --pheno ../../EUR/EUR.QC.pheno \
    --ignore-fid \
    --binary-target T \
    --thread 1 \
    --gtf ../../source_data/Homo_sapiens.GRCh37.87.gtf \
    --msigdb ../../source_data/h.all.v7.2.symbols.gmt \
    --multi-plot 10 \
    --out PRSet_BIP3_EUR_hg19_h.hallmark


# For Atopic Derm
#PRSet
# Msigdb C2: CGP (Chemical and genetic perturbations; n=3358 gene sets)
 Rscript ../../../PRSice_linux/PRSice.R --dir .\
    --prsice ../../../PRSice_linux/PRSice_linux \
    --base ../../../data/AtopicDerm.sumstats \
    --target ../../../data/PAPG-townsk05_1000gEUR_GRCh37.QC \
    --pheno ../../../data/pheno \
    --ignore-fid \
    --binary-target T \
    --thread 1 \
    --gtf ../../../source_data/Homo_sapiens.GRCh37.87.gtf \
    --msigdb ../../../source_data/c2.cgp.v7.2.symbols.gmt \
    --multi-plot 10 \
    --out AtopicDerm_PRSet_PAPG-townsk05_1000g_EUR_c2.cgp


# Msigdb C2: All curated gene sets (CGP & Canonical Pathways (CP) from BIOCARTA, KEGG, PID, REACTOME,and WikiPathways; n=6226 gene sets)
 Rscript ../../../PRSice_linux/PRSice.R --dir .\
    --prsice ../../../PRSice_linux/PRSice_linux \
    --base  ../../../data/AtopicDerm.sumstats \
    --target ../../../data/PAPG-townsk05_1000gEUR_GRCh37.QC \
    --pheno ../../../data/pheno \
    --ignore-fid \
    --binary-target T \
    --thread 1 \
    --gtf ../../../source_data/Homo_sapiens.GRCh37.87.gtf \
    --msigdb ../../../source_data/c2.all.v7.2.symbols.gmt \
    --multi-plot 10 \
    --out AtopicDerm_PRSet_PAPG-townsk05_1000g_EUR_c2.all

# Msigdb C5: Ontology gene sets ; n=14765 gene sets)
 Rscript ../../../PRSice_linux/PRSice.R --dir .\
    --prsice ../../../PRSice_linux/PRSice_linux \
    --base  ../../../data/AtopicDerm.sumstats \
    --target ../../../data/PAPG-townsk05_1000gEUR_GRCh37.QC \
    --pheno ../../../data/pheno \
    --ignore-fid \
    --binary-target T \
    --thread 1 \
    --gtf ../../../source_data/Homo_sapiens.GRCh37.87.gtf \
    --msigdb ../../../source_data/c5.all.v7.2.symbols.gmt \
    --multi-plot 10 \
    --out AtopicDerm_PRSet_PAPG-townsk05_1000g_EUR_c5.GO

# Msigdb C7: Immunologic signature gene sets; n=4872 gene sets)
 Rscript ../../../PRSice_linux/PRSice.R --dir .\
    --prsice ../../../PRSice_linux/PRSice_linux \
    --base  ../../../data/AtopicDerm.sumstats \
    --target ../../../data/PAPG-townsk05_1000gEUR_GRCh37.QC \
    --pheno ../../../data/pheno \
    --ignore-fid \
    --binary-target T \
    --thread 1 \
    --gtf ../../../source_data/Homo_sapiens.GRCh37.87.gtf \
    --msigdb ../../../source_data/c7.all.v7.2.symbols.gmt \
    --multi-plot 10 \
    --out AtopicDerm_PRSet_PAPG-townsk05_1000g_EUR_c7.immune

# Msigdb H: Hallmark gene sets; n=50 gene sets)
 Rscript ../../../PRSice_linux/PRSice.R --dir .\
    --prsice ../../../PRSice_linux/PRSice_linux \
    --base  ../../../data/AtopicDerm.sumstats \
    --target ../../../data/PAPG-townsk05_1000gEUR_GRCh37.QC \
    --pheno ../../../data/pheno \
    --ignore-fid \
    --binary-target T \
    --thread 1 \
    --gtf ../../../source_data/Homo_sapiens.GRCh37.87.gtf \
    --msigdb ../../../source_data/h.all.v7.2.symbols.gmt \
    --multi-plot 10 \
    --out AtopicDerm_PRSet_PAPG-townsk05_1000g_EUR_h.hallmark









