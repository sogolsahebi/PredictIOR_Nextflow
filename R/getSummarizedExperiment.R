#!/usr/bin/env Rscript

# Goal: 
#- Function to convert a MAE to multiple SummarizedExperiments based on all unique treatment types

# Load required packages
library(SummarizedExperiment)

getSummarizedExperiment <- function(maepath, outputPath) {
  # Check if the input file exists
  if (!file.exists(maepath)) {
    stop("Input file does not exist:", maepath)
  }

  # Read MAE from the specified path
  mae <- readRDS(maepath)
  
  # Extract components
  expr <- assays(mae)[["expr"]]
  clin <- data.frame(colData(mae))
  annot <- data.frame(rowData(mae@ExperimentList$expr))
  
  # Filter annotations for protein-coding genes
  annot <- annot[annot$gene_type == "protein_coding", ]
  annot <- annot[!duplicated(annot$gene_name), ]
  
  # Adjust row names
  rownames(expr) <- annot$gene_name
  
  # Process each unique treatment
  uniqueTreatments <- unique(clin$treatment)
  seList <- list()
  
  for (treatment in uniqueTreatments) {
    clinSubset <- clin[clin$treatment == treatment, ]
    exprSubset <- expr[, rownames(clinSubset), drop=FALSE]
    
    # Create a SummarizedExperiment
    se <- SummarizedExperiment(assays = list(gene_expression = exprSubset),
                               rowData = annot, colData = clinSubset)
    
    # File naming and saving
    fileName <- sprintf("%s/ICB_Ravi__%s__%s.rda", outputPath, clinSubset$cancer_type[1], gsub("[+/]", "_", treatment))
    if (!dir.exists(outputPath)) {
      dir.create(outputPath, recursive = TRUE)
    }
    save(se, file=fileName)
    seList[[treatment]] <- se
    message("Processed and saved data for treatment: ", treatment)
  }
  
  return(seList)
}

# Example usage
outputPath <- "results/"
maepath <- "data/ICB_Ravi.rds"
mae_SEs <- getSummarizedExperiment(maepath, outputPath)
