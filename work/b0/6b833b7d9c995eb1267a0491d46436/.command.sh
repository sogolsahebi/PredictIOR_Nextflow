#! /usr/bin/env Rscript
file_paths <- c("ICB_Ravi__Lung__PD-L1_expr.csv", "ICB_Ravi__Lung__PD-L1_clin.csv")
lapply(file_paths, function(f) {
  cat("Contents of file: ", f, "\n")
  print(read.csv(f, head = TRUE))  # Just print the header or first few lines for checking
})
