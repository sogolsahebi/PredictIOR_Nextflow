#!/usr/bin/env Rscript
source('/R/load_libraries.R')

expr <- read.csv("ICB_small_Mariathasan_expr.csv", row.names = 1)
clin <- as.data.frame(colData(ICB_small_Mariathasan))

cox_os <- geneSurvCont(
    dat.icb = 'ICB_small_Mariathasan',
    clin = NULL,
    time.censor = 36,
    missing.perc = 0.5,
    const.int = 0.001,
    n.cutoff = 15,
    feature = c("CXCL9", "CXCL10", "TIGIT", "CD83", "STAT1", "CXCL11", "CXCL13", "CD8A", "CTLA4"),
    study = "ICB_small_Mariathasan",
    surv.outcome = OS,
    cancer.type = "Bladder",
    treatment = "PD-1/PD-L1"
)

# Adjust p-values for multiple testing
cox_os$FDR <- p.adjust(cox_os$Pval, method = "BH")
cox_os <- cox_os[order(cox_os$FDR), ]
write.csv(cox_os, file = "ICB_small_Mariathasan_cox_os_genes.csv", row.names = FALSE)
