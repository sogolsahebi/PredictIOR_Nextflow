#!/usr/bin/env Rscript
library(SummarizedExperiment)

data_icb = load("ICB_small_Mariathasan.rda")

expr <- assay(dat_icb)
clin <- as.data.frame(colData(dat_icb))
annot <- as.data.frame(rowData(dat_icb))

write.csv(expr, "ICB_small_Mariathasan_expr.csv", row.names = TRUE)
write.csv(clin, "ICB_small_Mariathasan_clin.csv", row.names = FALSE)
write.csv(annot, "ICB_small_Mariathasan_annot.csv", row.names = TRUE)
