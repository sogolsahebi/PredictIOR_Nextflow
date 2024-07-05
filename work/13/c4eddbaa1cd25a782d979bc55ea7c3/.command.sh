#!/usr/bin/env Rscript
library(SummarizedExperiment)
library(PredictioR)

expr <- read.csv("ICB_small_Mariathasan_expr.csv", row.names = 1)
clin <- read.csv("ICB_small_Mariathasan_clin.csv")

genes <- PredictIO_Bareche$gene_name

geneSigScore <- geneSigGSVA(
                        dat.icb = expr,
                        clin = clin,
                        sig = CYT_Rooney,
                        sig.name = "CYT_Rooney",
                        missing.perc = 0.5,
                        const.int = 0.001,
                        n.cutoff = 15,
                        sig.perc = 0.8, 
                        study = "ICB_small_Mariathasan")
}
geneSigScore <- as.data.frame(geneSigScore)
# Save as CSV file
write.csv(geneSigScore <- geneSigScore, file = "ICB_small_Mariathasan_GeneSigScore.csv", row.names = FALSE)
