#!/usr/bin/env Rscript
library(SummarizedExperiment)
library(PredictioR)
library(ggplot2)

expr <- as.matrix(read.csv("ICB_small_Mariathasan_expr.csv", row.names = 1))
clin <- read.csv("ICB_small_Mariathasan_clin.csv", row.names = 1)

cancer_type <- names(table(clin$cancer_type)[table(clin$cancer_type) >= 15])

# Generate KM plot
KMPlot(
    status = clin$event_occurred_os[clin$cancer_type %in% cancer_type"],
    time = clin$survival_time_os[clin$cancer_type %in% cancer_type"],
    time.censor = 36,
    var = as.numeric(scale(expr["CXCL9", clin$cancer_type %in% "cancer_type"])),
    title = "CXCL9 and OS Association",
    xlab = "Time (Months)",
    ylab = "Overall Survival",
    n0.cutoff = 5,
    n1.cutoff = 5,
    var.type = FALSE
)

ggsave("ICB_small_Mariathasan_kmplot_os_CXCL9.png")
