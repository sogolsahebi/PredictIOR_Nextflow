#!/usr/bin/env Rscript
library(box)
box::use("/workspace/PredictIOR_Nextflow/R/getGeneAssociation.R")  # Ensure this is the correct absolute path

expr <- read.csv("ICB_Ravi__Lung__PD-L1_expr.csv", row.names = 1)
clin <- read.csv("ICB_Ravi__Lung__PD-L1_clin.csv")

# Extract parts of study ID for use in your analysis functions
study_id_parts <- unlist(strsplit("ICB_Ravi__Lung__PD-L1", "__"))
print(study_id_parts)

# Here you would typically call functions from your sourced R script
# Example function usage: result <- someFunctionFromYourScript(arguments)
# Remember to adapt this to your actual function names and arguments

# Save or process results as needed
