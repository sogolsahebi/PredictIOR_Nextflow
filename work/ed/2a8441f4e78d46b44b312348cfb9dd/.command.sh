#!/bin/bash -ue
Rscript -e "library(SummarizedExperiment);                 se <- readRDS('ICB_Ravi__Lung__PD-\(L\)1.rda');                 print(class(se))"
