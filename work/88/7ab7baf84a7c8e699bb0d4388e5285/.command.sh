#! /usr/bin/env Rscript
library(SummarizedExperiment)

load("ICB_Ravi__Lung__PD-L1.rda")

expr <- assay(dat_icb)
clin <- as.data.frame(colData(dat_icb))

output_dir <- "./results"
if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
}

expr_file_path <- paste0(output_dir, "/", "ICB_Ravi__Lung__PD-L1_expr.csv")
clin_file_path <- paste0(output_dir, "/", "ICB_Ravi__Lung__PD-L1_clin.csv")

write.csv(expr, expr_file_path, row.names = TRUE)
write.csv(clin, clin_file_path, row.names = FALSE)
