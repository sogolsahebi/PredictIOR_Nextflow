#!/bin/bash -ue
Rscript -e '
library(SummarizedExperiment)
se <- readRDS('ICB_Ravi__Lung__IO\\+combo.rda')
'
