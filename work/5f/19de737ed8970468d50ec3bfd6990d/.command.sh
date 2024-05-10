#!/bin/bash -ue
Rscript -e '
library(SummarizedExperiment)
if (!file.exists("ICB_Ravi__Lung__IO+combo.rda")) {
  stop("File not found: ICB_Ravi__Lung__IO+combo.rda")
}
try({
  se <- load("ICB_Ravi__Lung__IO+combo.rda")
  print(paste("Loaded objects:", paste(ls(), collapse=", ")))
}, silent=FALSE)
'
