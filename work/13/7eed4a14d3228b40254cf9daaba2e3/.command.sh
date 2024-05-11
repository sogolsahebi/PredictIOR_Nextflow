#!/bin/bash -ue
library(SummarizedExperiment)
library(dplyr)

# Load the .rda file, which directly loads into 'se'
load("ICB_Wolf__Breast__IO+chemo.rda")
se <- get(ls()[1])  # Get the first (or the only) object loaded, assumed to be a SummarizedExperiment

# Extract expression data and clinical data
expr <- assay(se)
clin <- as.data.frame(colData(se))

# Source the gene survival analysis script
source("/workspace/PredictIOR_Nextflow/R/getGeneAssociation.R")

# Apply logistic regression function over the first 100 genes in the dataset.
logreg <- lapply(1:100, function(i) {
    # Call the geneLogReg function for each gene with relevant parameters
    res <- geneLogReg(dat.icb = expr,
                      clin = clin,
                      missing.perc = 0.5, # Exclude genes with over 50% missing data
                      const.int = 0.001, # A small constant added for numerical stability
                      n.cutoff = 15, # Minimum sample count required
                      feature = rownames(expr)[i], # Current gene
                      study = paste(strsplit("ICB_Wolf__Breast__IO+chemo", "__")[[1]][1],
                                    strsplit("ICB_Wolf__Breast__IO+chemo", "__")[[1]][2],
                                    strsplit("ICB_Wolf__Breast__IO+chemo", "__")[[1]][3],
                                    sep="__"), # Study ID
                      n0.cutoff = 3, # Minimum samples in one group
                      n1.cutoff = 3, # Minimum samples in another group
                      cancer.type = strsplit("ICB_Wolf__Breast__IO+chemo", "__")[[1]][2], # Cancer type
                      treatment = strsplit("ICB_Wolf__Breast__IO+chemo", "__")[[1]][3]) # Treatment type

    # Filter results lacking a coefficient


    res
})


# Write results to CSV
write.csv(logreg, "logreg_results.csv", row.names = FALSE)
