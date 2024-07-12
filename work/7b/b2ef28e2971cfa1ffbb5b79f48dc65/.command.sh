#!/usr/bin/env Rscript
source('/R/load_libraries.R')

expr <- read.csv("ICB_small_Padron_expr.csv", row.names = 1)
clin <- read.csv("ICB_small_Padron_clin.csv")

logreg <- geneLogReg(
    dat.icb = expr,
    clin = clin,
    missing.perc = 0.5,
    const.int = 0.001,
    n.cutoff = 15,
    feature = c("CXCL9", "CXCL10", "TIGIT", "CD83", "STAT1", "CXCL11", "CXCL13", "CD8A", "CTLA4"),
    study = "ICB_small_Padron", 
    n0.cutoff = 3, 
    n1.cutoff = 3,
    cancer.type = "Pancreas",
    treatment = "PD-1/PD-L1"
)

# Adjust P-values and sort by FDR
logreg <- logreg[order(logreg$FDR <- p.adjust(logreg$Pval, method = "BH"))]

# Save as CSV file
write.csv(logreg, file = "ICB_small_Padron_Response.csv", row.names = FALSE)
