#!/usr/bin/env Rscript
expr <- read.csv("ICB_Ravi__Lung__PD-L1_expr.csv", row.names = 1)
clin <- read.csv("ICB_Ravi__Lung__PD-L1_clin.csv")

# Check if the source file can be accessed and the function can be sourced without error
source("/workspace/PredictIOR_Nextflow/bin/getGeneAssociation.R")

# Ensure the split and paste are working as expected
print([ICB_Ravi, Lung, PD-L1])

# Placeholder for actual function call
# Insert actual function call here once above checks are confirmed to be correct
