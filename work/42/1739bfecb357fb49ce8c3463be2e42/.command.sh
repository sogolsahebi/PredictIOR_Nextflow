#!/usr/bin/env Rscript
library(SummarizedExperiment)

load("ICB_small_Mariathasan.rda")

expr <- as.matrix(assay(ICB_small_Mariathasan))
clin <- as.data.frame(colData(ICB_small_Mariathasan))
annot <- as.data.frame(rowData(ICB_small_Mariathasan))

write.csv(expr, "ICB_small_Mariathasan_expr.csv", row.names = TRUE)
write.csv(clin, "ICB_small_Mariathasan_clin.csv", row.names = FALSE)
write.csv(annot, "ICB_small_Mariathasan_annot.csv", row.names = TRUE)
