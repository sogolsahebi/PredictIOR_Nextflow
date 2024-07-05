#!/usr/bin/env Rscript
library(SummarizedExperiment)
library(PredictioR)

expr <- read.csv("ICB_Ravi__Lung__PD-L1_expr.csv", row.names = 1)
clin <- read.csv("ICB_Ravi__Lung__PD-L1_clin.csv")

logreg <- lapply(1:100, function(i) {
  res <- geneLogReg(
    dat.icb = expr,
    clin = clin,
    missing.perc = 0.5,
    const.int = 0.001,
    n.cutoff = 15,
    feature = rownames(expr)[i],
    study = 'ICB_Ravi__Lung__PD-L1',
    n0.cutoff = 3,
    n1.cutoff = 3,
    cancer.type = 'Lung',
    treatment = 'PD-L1'
  )

  # Filter out the results that do not have a coefficient (Coef) value.
  res <- res[!is.na(res$Coef), ]
  res
})

# Convert list to a data frame
logreg_df <- do.call(rbind, logreg)

# Save the results to a CSV file
write.csv(logreg_df, file = "ICB_Ravi__Lung__PD-L1_logreg_results.csv", row.names = FALSE)