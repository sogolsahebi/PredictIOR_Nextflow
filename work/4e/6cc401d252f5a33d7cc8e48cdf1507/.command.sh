#!/usr/bin/env Rscript
library(SummarizedExperiment)
library(PredictioR)
library(GSVA)

# Load icb rda file
load("ICB_small_Mariathasan.rda")
expr <- assay(ICB_small_Mariathasan)
clin <- as.data.frame(colData(ICB_small_Mariathasan))

# Load signature genes
load("null.rda") #this gives us 'sig'

geneSigScore <- geneSigGSVA(dat.icb = expr,
                            clin = clin,
                            sig = sig,
                            sig.name = "CYT_Rooney",
                            missing.perc = 0.5,
                            const.int = 0.001,
                            n.cutoff = 15,
                            sig.perc = 0.8, 
                            study = "ICB_small_Mariathasan")

geneSigScore <- as.data.frame(geneSigScore)
# Save as CSV file
write.csv(geneSigScore, file = "ICB_small_Mariathasan_GeneSigScore_GSVA.csv", row.names = FALSE)
