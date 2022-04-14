#!/bin/bash

#SBATCH -J Regenie_binary_genomewide
#SBATCH -t 18:0:0
#SBATCH --mem=10G
#SBATCH --output=./errors/RegenieBinary%A_%a.log
#SBATCH -a 1-267


projectDir=


genotypeDir=
regenie=${projectDir}/regenie_v2.2.4.gz_x86_64_Centos7_mkl
genoDataDir=${projectDir}/geneticData
phenoDataDir=${projectDir}/phenotypes
resultsDir=${projectDir}/assocResults
tmp=${projectDir}/tempDir

cd $projectDir

ch=$SLURM_ARRAY_TASK_ID

read chr chunk start end < <(sed -n ${ch}p ${genoDataDir}/chunks_minMaf0.0001_minInfo0.8.txt)

rsync -av ${genotypeDir}/imputed/cleanedEuro_chr${chr}_chunk${chunk}.* ${genoDataDir}

$regenie \
    --step 2 \
    --bgen ${genoDataDir}/cleanedEuro_chr${chr}_chunk${chunk}.bgen \
    --sample ${genoDataDir}/cleanedEuro_chr${chr}_chunk${chunk}.sample \
    --ref-first \
    --test additive \
    --bsize 1000 \
    --phenoFile ${phenoDataDir}/phenotypesBinary.txt \
    --covarFile ${phenoDataDir}/covariates.txt \
    --bt \
    --firth \
    --firth-se \
    --approx \
    --use-null-firth ${resultsDir}/binaryTraits_step1_firth.list \
    --pred ${resultsDir}/binaryTraits_step1_pred.list \
    --out ${resultsDir}/chr${chr}_chunk${chunk}_binaryTraitsResults
