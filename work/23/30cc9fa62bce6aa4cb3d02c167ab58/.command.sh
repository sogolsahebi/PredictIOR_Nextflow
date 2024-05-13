#! /usr/bin/env Rscript
library(SummarizedExperiment)

# Load the .rda file
load("ICB_Ravi__Lung__PD-L1.rda")

# Extract expression data and clinical data
expr <- assay(dat_icb)
clin <- as.data.frame(colData(dat_icb))

# Set the output directory
output_dir <- "./results"

# Define file paths for expression and clinical data
expr_file_path <- file.path(output_dir, "ICB_Ravi__Lung__PD-L1_expr.csv")
clin_file_path <- file.path(output_dir, "ICB_Ravi__Lung__PD-L1_clin.csv")

# Save the expression and clinical data to CSV files
write.csv(expr, expr_file_path, row.names = TRUE)
write.csv(clin, clin_file_path, row.names = FALSE)

# Output the file paths for Nextflow
out_expr_csv = expr_file_path
out_clin_csv = clin_file_path
