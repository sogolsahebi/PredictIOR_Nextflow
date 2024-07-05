#!/usr/bin/env Rscript
library(SummarizedExperiment)
library(PredictioR)

expr <- read.csv("ICB_small_Mariathasan_expr.csv", row.names = 1)
clin <- read.csv("ICB_small_Mariathasan_clin.csv")
#TODO
# cancer_type <- names(table(clin$cancer_type)[table(clin$cancer_type) >= 15])

# load "PredictIO_Bareche.rda"
load("PredictIO_Bareche.rda")

genes <- PredictIO_Bareche$gene_name

logreg <- geneLogReg(dat.icb = expr,
                 missing.perc = 0.5,
                 const.int = 0.001,
                 n.cutoff = 15,
                 feature = genes,
                 study = "ICB_small_Mariathasan", 
                 n0.cutoff = 10,
                 n1.cutoff = 10,
                 cancer.type = cancer_type,
                 treatment = 'PD-1/PD-L1')

# Adjust P-values and sort by FDR
logreg <- logreg[order(logreg$FDR <- p.adjust(logreg$Pval, method = "BH")), ] 

# Save as CSV file
write.csv(logreg, file = "ICB_small_Mariathasan_Response_RvsNR.csv", row.names = FALSE)
