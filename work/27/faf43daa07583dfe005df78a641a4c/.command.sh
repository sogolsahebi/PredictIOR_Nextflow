#!/usr/bin/env Rscript
library(SummarizedExperiment)
library(PredictioR)

# Load icb rda file
load("ICB_small_Mariathasan.rda")
expr <- assay(ICB_small_Mariathasan)
clin <- as.data.frame(colData(ICB_small_Mariathasan))

# Load signature genes
load("EMT_Thompson.rda") // this gives us 'sig'

geneSigScore <- geneSigWeightedMean(dat.icb = expr,
                                    clin = clin,
                                    sig = sig,
                                    sig.name = "EMT_Thompson",
                                    missing.perc = 0.5,
                                    const.int = 0.001,
                                    n.cutoff = 15,
                                    sig.perc = 0.8, 
                                    study = "ICB_small_Mariathasan")

geneSigScore <- as.data.frame(geneSigScore)
# Save as CSV file
write.csv(geneSigScore, file = "ICB_small_Mariathasan_GeneSigScore_WeightedMean.csv", row.names = FALSE)
