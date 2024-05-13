#! /usr/bin/env Rscript

 library(SummarizedExperiment)

 if (!require("BiocManager", quietly = TRUE))
 install.packages("BiocManager")

 BiocManager::install("SummarizedExperiment")

 # Load the .rda file
 load("ICB_Ravi__Lung__PD-L1.rda")

 # Extract expression data and clinical data
 expr <- assay(se)
 clin <- as.data.frame(colData(se))

 # Define file paths dynamically based on parameters
 expr_file_path <- "./results/ICB_Ravi__Lung__PD-L1_expr.csv"
 clin_file_path <- "./results/ICB_Ravi__Lung__PD-L1_clin.csv"

 # Save the data frames for visibility and further analysis
 write.csv(expr, expr_file_path, row.names = TRUE)
 write.csv(clin, clin_file_path, row.names = FALSE)
