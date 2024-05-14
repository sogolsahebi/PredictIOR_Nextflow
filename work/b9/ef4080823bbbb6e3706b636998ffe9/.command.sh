#!/usr/bin/env Rscript
library(SummarizedExperiment)

args <- commandArgs(trailingOnly = TRUE)
study_id_parts <- unlist(strsplit(args[1], " "))

load("ICB_Ravi__Lung__PD-L1.rda")
expr <- assay(dat_icb)
clin <- as.data.frame(colData(dat_icb))

source('/R_scripts/getGeneAssociation.R')
source('/R_scripts/getHR.R')

results <- lapply(1:100, function(i) {
  geneSurvCont(
    dat.icb = expr,
    clin = clin,
    time.censor = 36,
    missing.perc = 0.5,
    const.int = 0.001,
    n.cutoff = 15,
    feature = rownames(expr)[i],
    study = paste(study_id_parts[1], study_id_parts[2], study_id_parts[3], sep='__'),
    surv.outcome = 'OS',
    cancer.type = study_id_parts[2],
    treatment = study_id_parts[3]
  )
})

write.csv(results, file = 'ICB_Ravi__Lung__PD-L1_results.csv', row.names = FALSE)
