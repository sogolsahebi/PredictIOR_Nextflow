#!/usr/bin/env Rscript
library(SummarizedExperiment)
library(PredictioR)
library(survival)

expr <- as.matrix(read.csv("ICB_small_Mariathasan_expr.csv", row.names = 1))
clin <- read.csv("ICB_small_Mariathasan_clin.csv", row.names = 1)

# Ensure cancer_type is extracted correctly
cancer_type <- names(table(clin$cancer_type)[table(clin$cancer_type) >= 15])

# Perform survival analysis
dicho <- survDicho( status = dat_clin$event_occurred_os ,
       time = dat_clin$survival_time_os,
       time.censor= 36,
       var = as.numeric(expr["CXCL9", ]),
       n0.cutoff = 5,
       n1.cutoff = 5,
       method = "median",
       var.type = FALSE)


# convert to dataframe
dicho <- as.data.frame(t(dicho))

write.csv(dicho, file = "ICB_small_Mariathasan_dicho_os_CXCL9.csv", row.names = FALSE)
