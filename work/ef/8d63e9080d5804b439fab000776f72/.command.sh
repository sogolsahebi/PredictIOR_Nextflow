#!/usr/bin/env Rscript
library(SummarizedExperiment)
library(PredictioR)
library(GSVA)

expr <- read.csv("ICB_small_Mariathasan_expr.csv", row.names = 1)
clin <- read.csv("ICB_small_Mariathasan_clin.csv")

geneSigScore_WeightedMean <- geneSigMean(dat.icb = expr,
                            clin = clin,
                            sig = EMT_Thompson,
                            sig.name = 'EMT_Thompson',
                            missing.perc = 0.5,
                            const.int = 0.001,
                            n.cutoff = 15,
                            sig.perc = 0.8, 
                            study = "ICB_small_Mariathasan")


geneSigScore_WeightedMean  <- as.data.frame(geneSigScore_WeightedMean)

# Save as CSV file
write.csv(geneSigScore_WeightedMean  file = "ICB_small_Mariathasan_GeneSigScore_WeightedMean.csv", row.names = FALSE)
