#!bin/bash
weight=$1
gwas=$2
prefix=$3
weight_dir=$4
FUSIONi=/mnt/sda/prog/fusion_twas-master
OUT=/mnt/sda/tsk/meta/finn_ukb_r9/panc/fusion/${weight}
mkdir $OUT

cd $OUT
mkdir -p $OUT/${prefix}

for chr in $(seq 12 12); do
Rscript $FUSIONi/FUSION.assoc_test.R  --sumstats $gwas --weights ${weight_dir}/${weight}.pos  --weights_dir ${weight_dir}/${weight}  --ref_ld_chr /mnt/sda/resource/ref/LD/1KG_202204/chr --chr ${chr} --out $OUT/${prefix}/chr${chr}.out --force_model ldpred2
done
