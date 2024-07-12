#!/usr/bin/env Rscript
source('/R/load_libraries.R')

load("ICB_small_Liu.rda")

expr <- assay(ICB_small_Liu)
clin <- as.data.frame(colData(ICB_small_Liu))
annot <- as.data.frame(rowData(ICB_small_Liu))

write.csv(expr, "ICB_small_Liu_expr.csv", row.names = TRUE)
write.csv(clin, "ICB_small_Liu_clin.csv", row.names = FALSE)
write.csv(annot, "ICB_small_Liu_annot.csv", row.names = TRUE)
