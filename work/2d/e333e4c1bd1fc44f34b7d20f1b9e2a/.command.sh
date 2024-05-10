#!/bin/bash -ue
Rscript -e 'library(SummarizedExperiment); load("ICB_Ravi__Lung__IO+combo.rda"); saveRDS(get(ls()[1]), file="ICB_Ravi__Lung__IO+combo.rds")'
