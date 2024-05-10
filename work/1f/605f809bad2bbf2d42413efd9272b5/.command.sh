#!/bin/bash -ue
Rscript -e 
'library(SummarizedExperiment); 
se <- load("ICB_Wolf__Breast__IO+chemo.rda")

'
