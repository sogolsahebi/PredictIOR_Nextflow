#!/usr/bin/env Rscript
source('/R/load_libraries.R')

load("ICB_small_Padron.rda")

expr <- assay(ICB_small_Padron)
clin <- as.data.frame(colData(ICB_small_Padron))
annot <- as.data.frame(rowData(ICB_small_Padron))

write.csv(expr, "ICB_small_Padron_expr.csv", row.names = TRUE)
write.csv(clin, "ICB_small_Padron_clin.csv", row.names = FALSE)
write.csv(annot, "ICB_small_Padron_annot.csv", row.names = TRUE)
