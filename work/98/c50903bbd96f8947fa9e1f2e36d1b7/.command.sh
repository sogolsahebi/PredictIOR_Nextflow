#!/usr/bin/env Rscript
library(SummarizedExperiment)
library(PredictioR)
library(GSVA)

print("Starting GeneSigScore_PredictIO process")

expr <- read.csv("ICB_small_Mariathasan_expr.csv", row.names = 1)
clin <- read.csv("ICB_small_Mariathasan_clin.csv")

print("Expression and clinical data loaded successfully.")

# load "PredictIO_Bareche.rda"
load("PredictIO_Bareche.rda")

print("PredictIO_Bareche data loaded successfully.")

geneSigScore_PredictIO <- geneSigPredictIO(dat.icb = expr,
                        clin = clin, 
                        sig = PredictIO_Bareche,
                        sig.name = 'PredictIO_Bareche',
                        missing.perc = 0.5,
                        const.int = 0.001,
                        n.cutoff = 15,
                        sig.perc = 0.8, 
                        study = "ICB_small_Mariathasan")

print("geneSigPredictIO function executed.")

# Check if the object is created successfully
print(geneSigScore_PredictIO)

# Convert to dataframe
geneSigScore_PredictIO <- as.data.frame(t(geneSigScore_PredictIO))

# Print first few rows of the dataframe
print(head(geneSigScore_PredictIO))

# Save as CSV file
write.csv(geneSigScore_PredictIO, file = "ICB_small_Mariathasan_GeneSigScore_PredictIO.csv", row.names = FALSE)

print("File written successfully.")
