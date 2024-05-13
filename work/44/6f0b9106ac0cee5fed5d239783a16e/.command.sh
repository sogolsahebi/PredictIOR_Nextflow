#!/usr/bin/env Rscript
library(SummarizedExperiment)
load("ICB_Ravi__Lung__PD-L1.rda")
expr <- assay(dat_icb)
clin <- as.data.frame(colData(dat_icb))
write.csv(expr, "ICB_Ravi__Lung__PD-L1_expr.csv", row.names = TRUE)
write.csv(clin, "ICB_Ravi__Lung__PD-L1_clin.csv", row.names = FALSE)
