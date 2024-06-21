#!/usr/bin/env Rscript
library(SummarizedExperiment)

load("ICB_Ravi__Lung__PD-L1.rda")

expr <- assay(dat_icb)
clin <- as.data.frame(colData(dat_icb))
annot <- as.data.frame(rowData(dat_icb))

write.csv(expr, "ICB_Ravi__Lung__PD-L1_expr.csv", row.names = TRUE)
write.csv(clin, "ICB_Ravi__Lung__PD-L1_clin.csv", row.names = FALSE)
write.csv(annot, "ICB_Ravi__Lung__PD-L1_annot.csv", row.names = TRUE)
