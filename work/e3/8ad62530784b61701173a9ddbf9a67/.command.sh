#!/usr/bin/env Rscript
library(SummarizedExperiment)
library(PredictioR)
library(GSVA)

# Load icb rda file
load("ICB_small_Mariathasan.rda")
expr <- assay(ICB_small_Mariathasan)
clin <- as.data.frame(colData(ICB_small_Mariathasan))

# Load signature genes
load("PredictIO_Bareche.rda") # this gives us 'sig'

geneSigScore_PredictIO <- geneSigPredictIO(dat.icb = expr,
                                           clin = clin, 
                                           sig = sig,
                                           sig.name = "PredictIO_Bareche",
                                           missing.perc = 0.5,
                                           const.int = 0.001,
                                           n.cutoff = 15,
                                           sig.perc = 0.8, 
                                           study = "ICB_small_Mariathasan")

# Convert to dataframe
geneSigScore_PredictIO <- as.data.frame(t(geneSigScore_PredictIO))

# Save as CSV file
write.csv(geneSigScore_PredictIO, file = "ICB_small_Mariathasan_GeneSigScore_PredictIO.csv", row.names = FALSE)
