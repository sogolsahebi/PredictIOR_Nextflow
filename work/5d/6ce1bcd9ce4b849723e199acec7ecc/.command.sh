#!/usr/bin/env Rscript
expr <- read.csv("ICB_Ravi__Lung__PD-L1_expr.csv", row.names = 1)
clin <- read.csv("ICB_Ravi__Lung__PD-L1_clin.csv")
# Ensure the R script path is correct; update the path according to your environment
source("/workspace/PredictIOR_Nextflow/R/getGeneAssociation.R")
study_id_parts <- unlist(strsplit("ICB_Ravi__Lung__PD-L1", "__"))
print(study_id_parts)
# Placeholder for actual function call, insert your function here
