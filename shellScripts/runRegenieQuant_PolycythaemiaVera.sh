#!/bin/bash

#SBATCH -J Regenie_quant_PolycythaemiaVera
#SBATCH -t 24:0:0
#SBATCH --mem=32G
#SBATCH --output=./errors/RegenieQuantPolycythaemiaVera.log


projectDir=


regenie=${projectDir}/regenie_v2.2.4.gz_x86_64_Centos7_mkl
genoDataDir=${projectDir}/geneticData
phenoDataDir=${projectDir}/phenotypes
resultsDir=${projectDir}/assocResults
tmp=${projectDir}/tempDir

cd $projectDir


plink \
  --bfile ${genoDataDir}/grmSNPsEuro \
  --keep ${phenoDataDir}/PolycythaemiaVera.txt \
  --mac 10 \
  --write-snplist \
  --out ${genoDataDir}/snps_PolycythaemiaVera



$regenie  \
  --step 1 \
  --bed ${genoDataDir}/grmSNPsEuro \
  --extract ${genoDataDir}/snps_PolycythaemiaVera.snplist \
  --phenoFile ${phenoDataDir}/phenotypesQuant.txt \
  --covarFile ${phenoDataDir}/covariates.txt \
  --keep ${phenoDataDir}/PolycythaemiaVera.txt \
  --bsize 1000 \
  --lowmem \
  --lowmem-prefix ${tmp}/quantTraits_tmp_PolycythaemiaVera \
  --out ${resultsDir}/quantTraits_step1_PolycythaemiaVera


nCh=$(wc -l chunks.txt | cut -d " " -f 1)

for ch in $(seq 1 $nCh)
do

  read gene chr start end start500kb end500kb chunk < <(sed -n ${ch}p chunks.txt)

  $regenie \
      --step 2 \
      --bgen ${genoDataDir}/cleanedEuro_chr${chr}_chunk${chunk}.bgen \
      --sample ${genoDataDir}/cleanedEuro_chr${chr}_chunk${chunk}.sample \
      --range ${chr}:${start500kb}-${end500kb} \
      --ref-first \
      --test additive \
      --bsize 1000 \
      --phenoFile ${phenoDataDir}/phenotypesQuant.txt \
      --covarFile ${phenoDataDir}/covariates.txt \
      --keep ${phenoDataDir}/PolycythaemiaVera.txt\
      --pred ${resultsDir}/quantTraits_step1_PolycythaemiaVera_pred.list \
      --out ${resultsDir}/${gene}_quantTraitsResults_PolycythaemiaVera

done
