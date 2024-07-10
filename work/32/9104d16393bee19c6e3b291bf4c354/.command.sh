#!/usr/bin/env Rscript
source('/R/load_libraries.R')

expr <- read.csv("ICB_small_Mariathasan_expr.csv", row.names = 1)
clin <- read.csv("ICB_small_Mariathasan_clin.csv")

cox_pfs <- geneSurvCont(
    dat.icb = expr,
    clin = clin,
    time.censor = 36,
    missing.perc = 0.5,
    const.int = 0.001,
    n.cutoff = 15,
    feature = c("CXCL9", "CXCL10", "TIGIT", "CD83", "STAT1", "CXCL11", "CXCL13", "CD8A", "CTLA4"),
    study = "ICB_small_Mariathasan",
    surv.outcome = 'PFS',
    cancer.type = "Bladder",
    treatment = "PD-1/PD-L1"
)

# Adjust p-values for multiple testing
cox_pfs$FDR <- p.adjust(cox_pfs$Pval, method = "BH")
cox_pfs <- cox_pfs[order(cox_pfs$FDR), ]
write.csv(cox_pfs, file = "ICB_small_Mariathasan_cox_pfs_genes.csv", row.names = FALSE)
