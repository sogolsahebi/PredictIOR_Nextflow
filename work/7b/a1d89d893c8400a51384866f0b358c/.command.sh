#!/usr/bin/env Rscript
library(SummarizedExperiment)
library(PredictioR)
library(survival)

expr <- as.matrix(read.csv("ICB_small_Mariathasan_expr.csv", row.names = 1))
clin <- read.csv("ICB_small_Mariathasan_clin.csv", row.names = 1)

# Ensure cancer_type is extracted correctly
cancer_type <- names(table(clin$cancer_type)[table(clin$cancer_type) >= 15])

# Perform survival analysis
dicho <- survCont(
    status = clin$event_occurred_os,
    time = clin$survival_time_os,
    time.censor = 36,
    var = as.numeric(expr[CXCL9, ])
)
# convert to dataframe
dicho <- as.data.frame(t(dicho))

write.csv(dicho, file = "ICB_small_Mariathasan_dicho_os_CXCL9.csv", row.names = FALSE)
