#!/usr/bin/env Rscript
library(SummarizedExperiment)

args &lt;- commandArgs(trailingOnly = TRUE)
study_id_parts &lt;- unlist(strsplit(args[1], " "))

expr &lt;- read.csv("ICB_Ravi__Lung__PD-L1_expr.csv", row.names = 1)
clin &lt;- read.csv("ICB_Ravi__Lung__PD-L1_clin.csv")

source('/PredictIOR_Nextflow/R/getGeneAssociation.R')
source('/PredictIOR_Nextflow/R/getHR.R')

logreg &lt;- lapply(1:min(200, nrow(expr)), function(i) {
  res &lt;- geneLogReg(
    dat.icb = expr,
    clin = clin,
    missing.perc = 0.5,
    const.int = 0.001,
    n.cutoff = 15,
    feature = rownames(expr)[i],
    study = paste(study_id_parts[0], study_id_parts[1], study_id_parts[2], sep='__'),
    n0.cutoff = 3,
    n1.cutoff = 3,
    cancer.type = study_id_parts[1],
    treatment = study_id_parts[2]
  )

  res
})
res

write.csv(res, file = "ICB_Ravi__Lung__PD-L1_results_logreg.csv", row.names = FALSE)
