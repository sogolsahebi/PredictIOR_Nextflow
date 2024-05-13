#!/usr/bin/env Rscript
expr <- read.csv("ICB_Ravi__Lung__PD-L1_expr.csv", row.names = 1)
clin <- read.csv("ICB_Ravi__Lung__PD-L1_clin.csv")

# Check if the source file can be accessed and the function can be sourced without error
box::use("R/getGeneAssociation.R[...]")


# Ensure the split and paste are working as expected
study_id_parts <- unlist(strsplit("ICB_Ravi__Lung__PD-L1", "__"))
print(study_id_parts)

# Placeholder for actual function call
# Insert actual function call here once above checks are confirmed to be correct
