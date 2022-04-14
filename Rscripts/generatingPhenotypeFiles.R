library(data.table)
library(magrittr)
library(tidyverse)
library(lubridate)

dataDir <- ""
phenoDir <- ""

data <- fread(paste0(dataDir,"ukb42082.tab"), select = c(1, 8836:8865,  8983:9312), quote="")


getPhenotype <- function(icd10, icd9=NA) {


  icd10cols <- names(data) %like% "f.40001|f.40002|f.41201|f.41202|f.41204." %>% which
  hasICD10 <- apply(data[,..icd10cols], 1, function(r) ifelse(any(r %like% icd10), 1, 0))

  if(!is.na(icd9)) {
    icd9cols    <- names(data) %like% "f.41203.|f.41205."  %>% which

    hasICD9 <- apply(data[,..icd9cols], 1, function(r) ifelse(any(r %like% icd9), 1, 0))

  }

  if(!is.na(icd9)) {
    hasPheno <- ifelse(hasICD10==1 | hasICD9==1, 1, 0)
  } else {
    hasPheno <- hasICD10
  }

  return(hasPheno)

}

PolycythaemiaVera <- getPhenotype(icd10="D45", icd9="2384")
EssentialThrombocythaemia <- getPhenotype(icd10="D473", icd9="23871")
PrimaryMyelofibrosis  <-  getPhenotype(icd10="D474", icd9="28983")

phenotypes <- cbind(data[,f.eid],
                    data[,f.eid],
                    PolycythaemiaVera,
                    EssentialThrombocythaemia,
                    PrimaryMyelofibrosis) %>%
  as.data.table %>%
  setnames(., c("V1", "V2"), c("FID", "IID"))


## write sample lists


PolycythaemiaVera <- phenotypes[PolycythaemiaVera==1, .(FID, IID)]
fwrite(PolycythaemiaVera, file = paste0(phenoDir,"PolycythaemiaVera.txt"), sep = "\t", na = "NA", quote = F)

EssentialThrombocythaemia <- phenotypes[EssentialThrombocythaemia==1, .(FID, IID)]
fwrite(EssentialThrombocythaemia, file = paste0(phenoDir,"EssentialThrombocythaemia.txt"), sep = "\t", na = "NA", quote = F)

PrimaryMyelofibrosis <- phenotypes[PrimaryMyelofibrosis==1, .(FID, IID)]
fwrite(PrimaryMyelofibrosis, file = paste0(phenoDir,"PrimaryMyelofibrosis.txt"), sep = "\t", na = "NA", quote = F)

healthy <- phenotypes[PrimaryMyelofibrosis==0 & EssentialThrombocythaemia==0 & PolycythaemiaVera==0, .(FID, IID)]
fwrite(healthy, file = paste0(phenoDir,"Healthy.txt"), sep = "\t", na = "NA", quote = F)

## if not in "healthy" list, and not a case, set phenotype to NA
phenotypesOut <- phenotypes[PolycythaemiaVera==0 & !(IID %in% healthy[,IID]), PolycythaemiaVera := NA] %>%
 .[EssentialThrombocythaemia==0 & !(IID %in% healthy[,IID]), EssentialThrombocythaemia := NA] %>%
  .[PrimaryMyelofibrosis==0 & !(IID %in% healthy[,IID]), PrimaryMyelofibrosis := NA]

cols <- names(phenotypesOut)[-c(1,2)]
phenotypesOut[, lapply(.SD, table), .SDcols=cols]

fwrite(phenotypesOut, file = paste0(phenoDir,"phenotypesBinary.txt"), sep = "\t", na = "NA", quote = F)

## read in cleaned haem traits and merge with ICD phenotypes
haem <- fread("/stornext/Bioinf/data/lab_bahlo/users/jackson.v/MLKL_inflamm/phenotypes/cleanedPhenotypes.txt") %>%
  select(., FID, IID, haemConc, haemocritPerc, erythrocyteCount, corpuscularVol, corpuscularHaem, corpuscularHaemConc)


# pheno <- phenotypes[haem, on = c("FID", "IID")]
# fwrite(pheno, file = paste0(phenoDir,"phenotypes.txt"))
