#!bin/bash
gene=$1
pqtl=$2
weight_dir=$3
OUT=$4
gwasN=$5
gtf="/mnt/sda/resource/ref/gtf/Homo_sapiens.GRCh38.100.gtf.h.geneANNO"
info="/mnt/sdb/resource/xQTL/eQTL/catalogue/eQTL-Catalogue_info.txt"
chr=$(awk -v OFS='\t' -F'\t' -v gene=$gene '{if($5==gene){print $1}}' $gtf)
mkdir $OUT/${gene}/
#eqtl2twas test
for i in `cut -f1 "/mnt/sda/tsk/TWMR/Panc/4expo.txt"`;do
ENSG=$(awk -v OFS='\t' -F'\t' -v gene=$gene '{if($5==gene){print $4}}' $gtf)
panelN=$(awk -v OFS='\t' -F'\t' -v panel=$i '{if($1==panel){print $7}}' $info)
if grep -w $ENSG ${weight_dir}/${i}.pos; then
  Rscript /mnt/sda/prog/eQTL_to_TWAS/FUSION.assoc_test.edit.R --sumstats $pqtl --weights ${weight_dir}/${i}.pos --weights_dir ${weight_dir}/${i} --ref_ld_chr /mnt/sda/resource/ref/LD/1KG_202204/chr --chr ${chr} --coloc_P 0.01 --out $OUT/${gene}/${i}.chr${chr}.out --GWASN $gwasN --PANELN $panelN
else
  echo "There is no $gene in ${i}."
fi
done