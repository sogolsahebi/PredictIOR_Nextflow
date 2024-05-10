#!/bin/bash -ue
Rscript -e 'library(SummarizedExperiment); se <- load("ICB_Ravi__Lung__IO+combo.rda"); saveRDS(se, file="ICB_Ravi__Lung__IO+combo.rda")'
