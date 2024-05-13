#! /usr/bin/env Rscript
library(SummarizedExperiment)

# Load the .rda file
load("ICB_Ravi__Lung__PD-L1.rda")

# Extract expression data and clinical data
expr <- assay(dat_icb)
clin <- as.data.frame(colData(dat_icb))

# Save the data frames for visibility and further analysis
write.csv(expr, "./results/ICB_Ravi__Lung__PD-L1_expr.csv", row.names = TRUE)
write.csv(clin, "./results/ICB_Ravi__Lung__PD-L1_clin.csv", row.names = FALSE)
