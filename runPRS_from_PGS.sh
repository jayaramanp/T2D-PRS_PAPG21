#!/bin/sh

ml plink2
ml plink


genome=$1
out_dir=$2
ancestry=$3

#create bed files from personal genome vcfs
echo "running plink2"
plink2 --vcf $genome.vcf --max-alleles 2 --rm-dup --make-bed --out $genome
#plink2 --vcf PAPG21-3.vcf.gz --max-alleles 2 --rm-dup --make-bed --out PAPG21-3_capstone

####change fam file to this:
echo "0"$'\t'"PAPG21-3_capstone"$'\t'"0"$'\t'"0"$'\t'"2"$'\t'"2"$'\n' > $genome".fam"

###rerun plink2
plink2 --bfile $genome --max-alleles 2 --rm-dup --make-bed --out $genome

echo "getting snp dups"
cut -f 2 $genome.bim | sort | uniq -d > snp.dups

echo "making bfile"
plink --bfile $genome --exclude snp.dups --make-bed --out $genome


#echo "calculate PRS score from PGS"
#calculate genome-wide PRS score for a BIG 5 Personality trait
#plink --bfile $genome --score PGS000036.txt 2 4 6 sum --out $out_dir"/"$genome"_PGS000036_T2D_"$ancestry".PRSscore"

#echo "here is the output!"
#awk '{print "PGS000036_T2D" "\t" $0}' $out_dir"/"$genome"_PGS000036_"$ancestry".PRSscore.profile" | sed '1d' > $out_dir"/"$genome"_PGS000036_"$ancestry".PRSscore"

