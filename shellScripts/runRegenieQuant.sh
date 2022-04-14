#!/bin/bash

#SBATCH -J Regenie_quant
#SBATCH -t 24:0:0
#SBATCH --mem=32G
#SBATCH --output=/vast/scratch/users/jackson.v/hepcidinRegGenes/errors/RegenieQuant.log


projectDir=/vast/scratch/users/jackson.v/hepcidinRegGenes


regenie=${projectDir}/regenie_v2.2.4.gz_x86_64_Centos7_mkl
genoDataDir=${projectDir}/geneticData
phenoDataDir=${projectDir}/phenotypes
resultsDir=${projectDir}/assocResults
tmp=${projectDir}/tempDir

cd $projectDir

$regenie  \
  --step 1 \
  --bed ${genoDataDir}/grmSNPsEuro \
  --phenoFile ${phenoDataDir}/phenotypesQuant.txt \
  --covarFile ${phenoDataDir}/covariates.txt \
  --bsize 1000 \
  --lowmem \
  --lowmem-prefix ${tmp}/quantTraits_tmp \
  --out ${resultsDir}/quantTraits_step1


nCh=$(wc -l chunks.txt | cut -d " " -f 1)

for ch in $(seq 1 $nCh)
do

  read gene chr start end start50kb end50kb chunk < <(sed -n ${ch}p chunks.txt)

  $regenie \
      --step 2 \
      --bgen ${genoDataDir}/cleanedEuro_chr${chr}_chunk${chunk}.bgen \
      --sample ${genoDataDir}/cleanedEuro_chr${chr}_chunk${chunk}.sample \
      --range ${chr}:${start50kb}-${end50kb} \
      --ref-first \
      --test additive \
      --bsize 1000 \
      --phenoFile ${phenoDataDir}/phenotypesQuant.txt \
      --covarFile ${phenoDataDir}/covariates.txt \
      --pred ${resultsDir}/quantTraits_step1_pred.list \
      --out ${resultsDir}/${gene}_quantTraitsResults

done
