#!/usr/bin/env Rscript
library(SummarizedExperiment)

load("ICB_Mariathasan__Bladder__PD-L1..rda")

expr <- assay(ICB_Mariathasan__Bladder__PD-L1.)
clin <- as.data.frame(colData(ICB_Mariathasan__Bladder__PD-L1.))
annot <- as.data.frame(rowData(ICB_Mariathasan__Bladder__PD-L1.))

write.csv(expr, "ICB_Mariathasan__Bladder__PD-L1._expr.csv", row.names = TRUE)
write.csv(clin, "ICB_Mariathasan__Bladder__PD-L1._clin.csv", row.names = FALSE)
write.csv(annot, "ICB_Mariathasan__Bladder__PD-L1._annot.csv", row.names = TRUE)
