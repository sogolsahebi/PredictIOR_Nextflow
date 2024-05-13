#! /usr/bin/env Rscript
library(SummarizedExperiment)

# Print the current working directory
print(paste("Current working directory:", getwd()))

# Load the .rda file
load("ICB_Ravi__Lung__PD-L1.rda")

# Extract expression data and clinical data
expr <- assay(dat_icb)
clin <- as.data.frame(colData(dat_icb))

# Set the output directory
output_dir <- "./results"

# Check if the directory exists, create it if it doesn't
if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
}

# Define file paths dynamically based on parameters
expr_file_path <- paste0(output_dir, "/", "ICB_Ravi__Lung__PD-L1_expr.csv")
clin_file_path <- paste0(output_dir, "/", "ICB_Ravi__Lung__PD-L1_clin.csv")

# Save the data frames for visibility and further analysis
write.csv(expr, expr_file_path, row.names = TRUE)
write.csv(clin, clin_file_path, row.names = FALSE)

# Print the first few lines of each file to the console for quick verification
print(paste("First few lines of", expr_file_path, ":"))
print(read.csv(expr_file_path, nrows = 5))
print(paste("First few lines of", clin_file_path, ":"))
print(read.csv(clin_file_path, nrows = 5))