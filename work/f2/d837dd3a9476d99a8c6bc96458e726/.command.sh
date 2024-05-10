#!/bin/bash -ue
Rscript -e '
library(SummarizedExperiment);
if (!file.exists("ICB_Ravi__Lung__IO+combo.rda")) {
  stop("File not found: ICB_Ravi__Lung__IO+combo.rda");
}
load("ICB_Ravi__Lung__IO+combo.rda");  // Ensure this is actually outputting a .rds file
se <- get(ls()[1]);  // Assume the first object is the SummarizedExperiment
saveRDS(se, file="ICB_Ravi__Lung__IO+combo.rds");
'
