#!/usr/bin/env Rscript
library(SummarizedExperiment)

args <- commandArgs(trailingOnly = TRUE)
study_id_parts <- unlist(strsplit(args[1], " "))

expr <- read.csv("ICB_Ravi__Lung__PD-L1_expr.csv", row.names = 1)
clin <- read.csv("ICB_Ravi__Lung__PD-L1_clin.csv")

source('/R_scripts/getGeneAssociation.R')
source('/R_scripts/getHR.R')

# TODO: gene_level , adding Signiture level?
# for single gene?

results <- lapply(1:100, function(i) {
  geneSurvCont(
    dat.icb = expr,
    clin = clin,
    time.censor = 36,
    missing.perc = 0.5,
    const.int = 0.001,
    n.cutoff = 15,
    feature = rownames(expr)[i],
    study = 'ICB_Ravi__Lung__PD-L1',
    surv.outcome = 'OS',
    cancer.type = 'Lung',
    treatment = 'PD-L1'
  )
})

# Convert list to a data frame
results_df <- do.call(rbind, results)

write.csv(results_df, file = "ICB_Ravi__Lung__PD-L1_cox_os_results.csv", row.names = FALSE)
