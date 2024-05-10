#!/bin/bash -ue
Rscript -e 
'library(SummarizedExperiment); 
se <- load("ICB_Wolf__Breast__IO+chemo.rda")
expr <- assay(dat_icb)
clin <- colData(dat_icb) %>% as.data.frame()
'
