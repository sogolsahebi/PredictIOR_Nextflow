#!/usr/bin/env Rscript
library(SummarizedExperiment)

args <- commandArgs(trailingOnly = TRUE)
study_id_parts <- unlist(strsplit(args[1], "__"))

expr <- read.csv("ICB_Ravi__Lung__PD-L1_expr.csv", row.names = 1)
clin <- read.csv("ICB_Ravi__Lung__PD-L1_clin.csv")

source('/R_scripts/getGeneAssociation.R')
source('/R_scripts/getHR.R')

logreg <- lapply(1:min(200, nrow(expr)), function(i) {
  res <- geneLogReg(
    dat.icb = expr,
    clin = clin,
    missing.perc = 0.5,
    const.int = 0.001,
    n.cutoff = 15,
    feature = rownames(expr)[i],
    study = paste(study_id_parts[0], study_id_parts[1], study_id_parts[2], sep = '__'),
    n0.cutoff = 3,
    n1.cutoff = 3,
    cancer.type = study_id_parts[1],
    treatment = study_id_parts[2]
  )
  res
})

# Convert list to a data frame
logreg_df <- do.call(rbind, logreg)

# Save the results to a CSV file
write.csv(logreg_df, file = "ICB_Ravi__Lung__PD-L1_results_logreg.csv", row.names = FALSE)
