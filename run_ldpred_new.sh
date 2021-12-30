#########
### New Location
### /sc/arion/work/parani01

#1
#### Copy files from local to minerva
#scp /Users/iparanjpe/Documents/mount_sinai_research/ipm/heart_failure_prs/ldpred_input_files/ischemic_hf/ldpred_input_ischemic.ssv parani01@chimera.hpc.mssm.edu:/sc/orga/work/parani01/hf_prs/ischemic_hf/ldpred/input_files
# cp -r /sc/orga/work/parani01/kidney_stone_prs /hpc/users/parani01/parani01/

module load ldpred


base="/sc/arion/projects/bsr2401/PAPG21_22/jayarp02/capstone/"

target_file="/sc/private/regen/data/Regeneron/GSA/imputed_tgp_p3_plink/GSA_chr_all"

ref_file="T2D_TranEthnic.BMIadjusted.liftoverhg38.txt"
ldpred="/hpc/packages/minerva-centos7/ldpred/1.0.6/ldpred/LDpred.py"


summary_stats=${base}/ldpred/input_files/ldpred_input_ischemic.ssv
coord_output_file=${base}/ldpred/output_files/ldlc_ldpred1_standard.hdf5
ldpred_output_file=${base}/ldpred/output_files/ldlc_ldpred2_standard_out
ldlc_ldpred2_standard_ldfile=${base}/ldpred/output_files/local_ld

score_output_file="${base}/ldpred/output_files/score"


#1

nohup python $ldpred coord --gf=$ref_file --ssf=$summary_stats --N=500000  --out=$coord_output_file --ssf-format=STANDARD  &


#2
nohup python $ldpred gibbs --cf=$coord_output_file --ldr=257 --ldf=$ldlc_ldpred2_standard_ldfile --N=500000 --out=$ldpred_output_file >gibbs.log  &


#Job <125820421> is submitted to queue <private>.

#2200 sec

#R=120250/3000=round(40.08333)=40

#N=82974


####### Find duplicates to exclude
# cut -d " " -f 2 $target_file.bim | sort | uniq -d > duplicate_snps



duplicate_snps="/sc/orga/work/parani01/hf_prs/ischemic_hf/ldpred/duplicate_snps"

### run plink score function
module load plink

coef_file=${base}/ldpred/output_files/ldlc_ldpred2_standard_out_LDpred-inf.txt

nohup plink --bfile $target_file --exclude $duplicate_snps --score $coef_file 3 4 7 header  --out _inf &
