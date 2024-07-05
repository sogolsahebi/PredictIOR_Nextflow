#!/usr/bin/env Rscript
library(SummarizedExperiment)
library(PredictioR)
library(GSVA)

expr <- read.csv("ICB_small_Mariathasan_expr.csv", row.names = 1)
clin <- read.csv("ICB_small_Mariathasan_clin.csv")

# load "PredictIO_Bareche.rda"
load("PredictIO_Bareche.rda")

geneSigScore_PredictIO <- geneSigPredictIO(dat.icb = expr,
                            clin = clin, 
                            sig = PredictIO_Bareche,
                            sig.name = 'PredictIO_Bareche',
                            missing.perc = 0.5,
                            const.int = 0.001,
                            n.cutoff = 15,
                            sig.perc = 0.8, 
                            study = "ICB_small_Mariathasan")


geneSigScore_PredictIO  <- as.data.frame(t(geneSigScore_PredictIO))

# Save as CSV file
write.csv(geneSigScore_PredictIO, file = "ICB_small_Mariathasan_geneSigScore_PredictIO.csv", row.names = FALSE)
