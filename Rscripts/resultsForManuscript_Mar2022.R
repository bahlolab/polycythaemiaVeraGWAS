library(data.table)
library(magrittr)
library(tidyverse)
library(knitr)
library(qqman)
library(patchwork)
library(ggrepel)
library(cowplot)
library(boot)

source("myFunctions.R")

## gwas results for PV - additive model
chunkInfo <- fread("chunks_minMaf0.0001_minInfo0.8.txt", col.names = c( "chr", "chunk", "start", "end")) %>%
  .[, chunkID := paste0("chr",chr,"_chunk",chunk)]

chunks <- chunkInfo %>%
  .[c(1:267), chunkID]

tr <- "PolycythaemiaVera"
nCases <- 440

AFmin <- 10/(nCases*2)
AFmax <- 1-(10/(nCases*2))


results <- lapply(chunks, function(chunk) {

  prefix <- paste0("./assocResults/",chunk,"_binaryTraitsResults_PolycythaemiaVera")


  trResult <- paste0(prefix,".regenie") %>%
    fread(., check.names=T) %>%
    .[, P := 10^-LOG10P] %>%
    .[P==0, P := .Machine$double.xmin] %>%
    .[, c("N", "TEST", "EXTRA") := NULL] %>%
    setkey(., "CHROM", "GENPOS", "ID", "ALLELE0", "ALLELE1") %>%
    .[A1FREQ %between% c(AFmin, AFmax)]


  return(trResult)

}) %>%
  rbindlist



## recessive results - rs1800562 only, for each trait


traits <- c("haemConc",
            "haemocritPerc",
            "erythrocyteCount",
            "corpuscularVol",
            "corpuscularHaem",
            "corpuscularHaemConc")

analyses <- c("Healthy",
              "PolycythaemiaVera")



recessiveResults <- lapply(traits, function(tr) {

  stratResults <- lapply(analyses, function(set) {

     prefix <- paste0("./assocResults/HFE_quantTraitsRecessiveResults_",set,"_",tr)

    paste0(prefix,".regenie") %>%
      fread(., check.names=T) %>%
      .[ID == "rs1800562"] %>%
      .[, trait := tr] %>%
      .[, P := 10^-LOG10P] %>%
      .[P==0, P := .Machine$double.xmin] %>%
      .[, .(trait, ID, ALLELE0, ALLELE1, A1FREQ, BETA, SE, LOG10P, P)] %>%
      setnames(., c("A1FREQ", "BETA", "SE", "LOG10P", "P"), c(paste0("A1FREQ_",set), paste0("BETA_",set), paste0("SE_",set), paste0("LOG10P_",set), paste0("P_",set))) %>%
      setkey(., "trait", "ID", "ALLELE0", "ALLELE1")

  }) %>%
    Reduce(merge, .)

}) %>%
  rbindlist


fwrite(recessiveResults, file = "rs1800562_haemTraitsRecessive.csv", sep="\t")
