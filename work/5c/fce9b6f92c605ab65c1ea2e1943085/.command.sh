#!/bin/bash -ue
Rscript -e '
library(SummarizedExperiment);
if (!file.exists("ICB_Ravi__Lung__IO+combo.rds")) {
  stop("File not found: ICB_Ravi__Lung__IO+combo.rds")
}
se <- readRDS("ICB_Ravi__Lung__IO+combo.rds");  // Changed from load() to readRDS()
saveRDS(se, file="ICB_Ravi__Lung__IO+combo.rds");
'
