#!/usr/bin/env Rscript
source('/R/load_libraries.R')

rda_file <- "./ICB_data/ICB_small_Mariathasan.rda"
load(rda_file)

# Assuming the loaded R object's name corresponds to the study_id variable
study <- get(study_id)

expr <- assay(study)
clin <- as.data.frame(colData(study))
annot <- as.data.frame(rowData(study))

write.csv(expr, "ICB_small_Mariathasan_expr.csv", row.names = TRUE)
write.csv(clin, "ICB_small_Mariathasan_clin.csv", row.names = FALSE)
write.csv(annot, "ICB_small_Mariathasan_annot.csv", row.names = TRUE)
