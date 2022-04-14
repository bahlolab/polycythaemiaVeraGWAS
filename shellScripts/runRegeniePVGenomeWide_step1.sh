#!/bin/bash

#SBATCH -J Regenie_binary
#SBATCH -t 24:0:0
#SBATCH --mem=32G
#SBATCH --output=./errors/RegenieBinary.log


projectDir=/


regenie=${projectDir}/regenie_v2.2.4.gz_x86_64_Centos7_mkl
genoDataDir=${projectDir}/geneticData
phenoDataDir=${projectDir}/phenotypes
resultsDir=${projectDir}/assocResults
tmp=${projectDir}/tempDir

cd $projectDir

$regenie  \
  --step 1 \
  --bed ${genoDataDir}/grmSNPsEuro \
  --phenoFile ${phenoDataDir}/phenotypesBinary.txt \
  --covarFile ${phenoDataDir}/covariates.txt \
  --bt \
  --write-null-firth \
  --bsize 1000 \
  --lowmem \
  --lowmem-prefix ${tmp}/binaryTraits_tmp \
  --out ${resultsDir}/binaryTraits_step1
