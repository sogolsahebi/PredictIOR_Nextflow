#!/usr/bin/env Rscript
library(SummarizedExperiment)
library(PredictioR)

expr <- as.matrix(read.csv("ICB_small_Mariathasan_expr.csv", row.names = 1))
clin <- read.csv("ICB_small_Mariathasan_clin.csv")

# Ensure cancer_type is extracted correctly
cancer_type <- names(table(clin$cancer_type)[table(clin$cancer_type) >= 15])

# Perform survival analysis
cox_os <- survCont(
    status = clin$event_occurred_os,
    time = clin$survival_time_os,
    time.censor = 36,
    var = as.numeric(scale(expr["CXCL9", clin$sample_id]))
)

# convert to dataframe
cox_os <- as.data.frame(t(cox_os))

write.csv(cox_os, file = "ICB_small_Mariathasan_cox_os_CXCL9.csv", row.names = FALSE)