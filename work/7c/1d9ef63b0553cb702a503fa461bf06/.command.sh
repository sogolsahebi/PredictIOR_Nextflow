#!/bin/bash -ue
Rscript -e '
library(SummarizedExperiment);
if (!file.exists("ICB_Ravi__Lung__IO+combo.rda")) {
  stop("File not found: ICB_Ravi__Lung__IO+combo.rda");
}
# Assuming the loaded RDA file directly contains the SummarizedExperiment object
load("ICB_Ravi__Lung__IO+combo.rda");
se <- get(ls()[1]);  // Assuming the first object is the SummarizedExperiment
saveRDS(se, file="ICB_Ravi__Lung__IO+combo.rds");
'
