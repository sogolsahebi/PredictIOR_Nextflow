#! /usr/bin/env Rscript
library(data.table)

expr_data <- fread("ICB_Ravi__Lung__PD-L1_expr.csv")
clin_data <- fread("ICB_Ravi__Lung__PD-L1_clin.csv")

print("Expression Data (first few rows):")
print(head(expr_data))

print("Clinical Data (first few rows):")
print(head(clin_data))
