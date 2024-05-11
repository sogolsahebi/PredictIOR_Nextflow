#!/usr/bin/env nextflow
nextflow.enable.dsl=2

process ProcessGeneLogReg {
    container 'nextflow-r-env:latest'

    input:
    path rda_file

    output:
    path "logreg_results.csv"

    script:
    """
    library(SummarizedExperiment)
    library(dplyr)

    # Load the .rda file, which directly loads into 'se'
    load("${rda_file}")
    se <- get(ls()[1])  

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
                          study = paste(strsplit("ICB_Ravi__Lung__IO+combo", "__")[[1]][1],
                                        strsplit("ICB_Ravi__Lung__IO+combo", "__")[[1]][2],
                                        strsplit("ICB_Ravi__Lung__IO+combo", "__")[[1]][3],
                                        sep="__"), # Study ID
                          n0.cutoff = 3, 
                          n1.cutoff = 3, 
                          cancer.type = strsplit("ICB_Ravi__Lung__IO+combo", "__")[[1]][2], # Cancer type
                          treatment = strsplit("ICB_Ravi__Lung__IO+combo", "__")[[1]][3]) # Treatment type


        # Filter results lacking a coefficient
        print(class(res))
        print(head(res))

        res <- res[!is.na(res$Coef), ]
        res
    })

    # Combine results into a single data frame
    logreg <- do.call(rbind, logreg)

    # Filter out rows without gene names
    logreg <- logreg[!is.na(logreg$Gene), ]

    # Adjust p-values for multiple testing using the Benjamini-Hochberg method
    logreg$FDR <- p.adjust(logreg$Pval, method = "BH")

    # Write results to CSV
    write.csv(logreg, "logreg_results.csv", row.names = FALSE)
    """
}

workflow {
    rdaFile = file("${baseDir}/data/ICB_Ravi__Lung__IO+combo.rda")
    ProcessGeneLogReg(rdaFile)
}
