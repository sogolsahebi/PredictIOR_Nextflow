#!/usr/bin/env Rscript
source('/R/load_libraries.R')

expr <- read.csv("ICB_small_Mariathasan_expr.csv", row.names = 1)
clin <- read.csv("ICB_small_Mariathasan_clin.csv")

logreg <- geneLogReg(
    dat.icb = expr,
    clin = clin,
    missing.perc = 0.5,
    const.int = 0.001,
    n.cutoff = 15,
    feature = genes,
    study = genes, 
    n0.cutoff = 5, 
    n1.cutoff = 5,
    cancer.type = "Bladder",
    treatment = "PD-1/PD-L1"
)

# Adjust P-values and sort by FDR
logreg <- logreg[order(logreg$FDR <- p.adjust(logreg$Pval, method = "BH"))]

# Save as CSV file
write.csv(logreg, file = "ICB_small_Mariathasan_Response.csv", row.names = FALSE)
