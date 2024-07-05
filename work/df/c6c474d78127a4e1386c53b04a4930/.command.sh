#!/usr/bin/env Rscript
library(SummarizedExperiment)
library(PredictioR)
library(GSVA)

expr <- read.csv("ICB_small_Mariathasan_expr.csv", row.names = 1)
clin <- read.csv("ICB_small_Mariathasan_clin.csv")

# Load signature genes
sig <- CYT_Rooney

geneSigScore <- geneSigGSVA(dat.icb = expr,
                        sig = sig,
                        sig.name = "CYT_Rooney",
                        missing.perc = 0.5,
                        const.int = 0.001,
                        n.cutoff = 15,
                        sig.perc = 0.8, 
                        study = "ICB_small_Mariathasan")

geneSigScore <- as.data.frame(geneSigScore)
# Save as CSV file
write.csv(geneSigScore, file = "ICB_small_Mariathasan_GeneSigScore.csv", row.names = FALSE)
