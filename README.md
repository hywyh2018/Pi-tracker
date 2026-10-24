# Pi-tracker
Protein-identity tracker (Pi-tracker/π-tracker) is a framework to trace the source and target characteristics of plasma proteins.
## Table of contents
- [Introduction](#Introduction)
- [Installation](#Installation)
- [Quick tutorial](#Tutorial)
- [Input files](#Inputs)
- [Output files](#Outputs)
- [Contact](#Contact)
- [Notes](#Notes)
## Introduction
Investigations of plasma proteins are emerging as a promising strategy for identifying biomarkers and druggable targets. However, the challenges remain in tracking the source and molecular target of plasma proteins. Unravelling their tissue/cell-specific manifestations is critical for understanding of disease etiology and guiding therapeutic development. Pi-tracker is a framework to trace the source and target characteristics of plasma proteins, which utilized MR and PWAS by inferring the causal relationships between plasma proteins and other endophenotypes. It is developed and maintained by Yahui Wang at Zhuan Liao's Group.

In this tool, we

- integrate MR and PWAS tools to trace the source and target characteristics of plasma proteins
- integrate eQTL_to_TWAS package to train eQTL weight models
- enable parallel computation to train eQTL weight models / perform PWAS and colocalization analysis simultaneously for multiple genes.

## Installation
We recommend use mamba for installation:
```
mamba create -n Pi-tracker r-coloc r::r-mad=0.8_2.1 dnachun::r-plink2r r-base r-optparse=1.7.3 r-data.table=1.14.8 r-susieR=0.12.35 r-bigsnpr=1.12.2 r-remotes=2.4.2.1 r-dplyr=1.1.3 r-fdrtool=1.2.17 r-rcolorbrewer=1.1_3 r-r.utils=2.12.2 r-circlize bioconductor-complexheatmap
```
## Quick tutorial
### Example
We identified the origin of 9 candidate plasma proteins using tissue eQTL from GTEx or TwinUK, and cell-type eQTL from OneK1K or Perez2022. These eQTL datasets were used as the exposures to perform MR by the GSMR method, and used as the input traits to build the PWAS prediction models by the eQTL_to_TWAS method. The pQTL data of plasma proteins from deCODE, UKB and OMICSCIENCE were used as the outcomes in MR or target traits in PWAS. Finally, we calculate an omnibus test P value using the aggregated Cauchy association test (ACAT-O) with MR and PWAS P values, which could possibly improve the power while still managing to control for type I error.

MR test:
```
/mnt/sda/prog/gcta/gcta-1.94.1 --mbfile 1.ref.txt --gsmr-file 2.expo.txt 3.outcome.txt --gsmr-direction 0 --out GATA3.PWMR.R005.P001 --effect-plot --gsmr-snp-min 2 --heidi-thresh 0.05 --gwas-thresh 0.01 --clump-r2 0.05 --thread-num 10 --diff-freq 1 --maf 0.005
```

PWAS test:
```
nohup bash "/mnt/sda/tsk/TWMR/pi-tracker-pwas.sh" "CLEC4C" "/mnt/sda/resource/deCODE/Olink/pwas/UKBPPP.CLEC4C.fusion.txt" "/mnt/sda/resource/TWAS/eQTL_to_TWAS/eQTL-Catalogue/" "/mnt/sda/tsk/TWMR/Panc/" 46663 > /mnt/sda/tsk/TWMR/Panc/log/UKBPPP.PWAS.CLEC4C.log &
```
## Input files
Input includes:

- predicted cause sumstatFile: summary statistics files of xQTL
- predicted effect sumstatFile: summary statistics files of pQTL
- sumstatFile: summary statistics files of endophenotypes
- gtfFile: gtf annotation file
- refFile: reference panel file
- weights_dir: weight model directory
- weights: weight model file

## Output files
The pipeline will return three files:

- gene_name.cell_name.gsmr: MR results
- gene_name.cell_name.out: PWAS results
- gene_name.cell_name.txt: ACAT results

The GSMR output includes:

- gene: gene id
- snp: snp id
- weight: corresponding non-zero weight

The PWAS output includes:

- gene: gene id
- window: the chosen 3D window
- phi: the chosen penalty factor for epigenetic annotated variants
- r2: r2 estimated by summary-statistics based cross-validation

The ACAT-O output includes:

- gene: gene id
- snp: snp id
- weight: corresponding non-zero weight
## Contact
For any issues regarding Pi-tracker, please contact wangyahui@ncpsb.org.cn.
## Notes
Pi-tracker is based on GCTA, FUSION, eQTL_to_TWAS, and ACAT-O. We really appreciate the authors of these tools.
