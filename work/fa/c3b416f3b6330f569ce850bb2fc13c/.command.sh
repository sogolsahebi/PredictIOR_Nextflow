#!/usr/bin/env Rscript
    library(SummarizedExperiment)
    library(PredictioR)

    expr <- read.csv("ICB_small_Mariathasan_expr.csv", row.names = 1)
    clin <- read.csv("ICB_small_Mariathasan_clin.csv")

    cox_os <- geneSurvCont(
        dat.icb = ICB_small_Mariathasan,
        time.censor = 36,
        missing.perc = 0.5,
        const.int = 0.001,
        n.cutoff = 15,
        feature = "CXCL9",
        study = 'ICB_Mariathasan',
        surv.outcome = 'OS',
        cancer.type = 'Bladder',
        treatment = 'PD-1/PD-L1'
)

    #cox_os <- do.call(rbind, cox_os)

    # Additional filtering 
    # cox_os <- cox_os[!is.na(cox_os$Gene), ]
    # cox_os$FDR <- p.adjust(cox_os$Pval, method = "BH")

    write.csv(cox_os, file = "ICB_small_Mariathasan_cox_os_genes.csv", row.names = FALSE)
