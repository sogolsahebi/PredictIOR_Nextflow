#!/usr/bin/env Rscript
library(SummarizedExperiment)
library(PredictioR)

expr <- read.csv("ICB_small_Mariathasan_expr.csv", row.names = 1)
clin <- read.csv("ICB_small_Mariathasan_clin.csv")

# Ensure cancer_type is extracted correctly
cancer_type <- names(table(clin$cancer_type)[table(clin$cancer_type) >= 15])

# Filter data to include only rows with non-missing values for the relevant columns
clin_filtered <- clin[!is.na(clin$event_occurred_os) & !is.na(clin$survival_time_os) & clin$cancer_type %in% cancer_type, ]
expr_filtered <- expr[, colnames(expr) %in% clin_filtered$sample_id]
cancer_type <- names( table( clin$cancer_type )[ table(clin$cancer_type ) >= 15 ] )

# Perform survival analysis
cox_os <- survCont(
    status = clin_filtered$event_occurred_os,
    time = clin_filtered$survival_time_os,
    time.censor = 36,
    var = as.numeric(scale(expr_filtered["CXCL9", clin_filtered$sample_id]))
)

write.csv(cox_os, file = "ICB_small_Mariathasan_cox_os_with_CXCL9.csv", row.names = FALSE)
