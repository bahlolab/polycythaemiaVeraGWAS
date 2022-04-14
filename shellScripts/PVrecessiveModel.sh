#!/bin/bash

#SBATCH -J Regenie_binary
#SBATCH -t 24:0:0
#SBATCH --mem=32G
#SBATCH --output=./errors/RegenieBinary.log


projectDir=


regenie=${projectDir}/regenie_v2.2.4.gz_x86_64_Centos7_mkl
genoDataDir=${projectDir}/geneticData
phenoDataDir=${projectDir}/phenotypes
resultsDir=${projectDir}/assocResults
tmp=${projectDir}/tempDir

cd $projectDir

# step 1 already run, for additive model

ch=2

  read gene chr start end start500kb end500kb chunk < <(sed -n ${ch}p chunks.txt)

  $regenie \
      --step 2 \
      --bgen ${genoDataDir}/cleanedEuro_chr${chr}_chunk${chunk}.bgen \
      --sample ${genoDataDir}/cleanedEuro_chr${chr}_chunk${chunk}.sample \
      --range ${chr}:${start500kb}-${end500kb} \
      --ref-first \
      --test recessive \
      --bsize 1000 \
      --phenoFile ${phenoDataDir}/phenotypesBinary.txt \
      --phenoCol PolycythaemiaVera \
      --covarFile ${phenoDataDir}/covariates.txt \
      --bt \
      --firth \
      --firth-se \
      --approx \
      --use-null-firth ${resultsDir}/binaryTraits_step1_firth.list \
      --pred ${resultsDir}/binaryTraits_step1_pred.list \
      --out ${resultsDir}/${gene}_PVRecessiveResults
