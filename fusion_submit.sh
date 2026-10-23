#!/bin/bash
prefix=$1
gwas=$2
weight_dir=$3
mkdir /mnt/sda/tsk/meta/secBM/fusion/${prefix}
cd /mnt/sda/resource/TWAS/FUSION/${weight_dir}/
for i in `ls *.pos`;do
weight=$(echo ${i//\.pos/})
bash /mnt/sda/tsk/meta/secBM/fusion/fusion_run.sh $weight $gwas $prefix $weight_dir
done
