#!/usr/bin/env nextflow
nextflow.enable.dsl=2

params.study_id = 'ICB_Ravi__Lung__PD-(L)1'
params.data_dir = './data'
params.out_dir = './results'
params.rda_file = "${params.data_dir}/${params.study_id}.rda"

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
    sedat_icb <- load("${params.rda_file}")  

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




#####################################################################

#!/usr/bin/env nextflow

params.study_id = 'ICB_Ravi__Lung__PD-(L)1'
params.data_dir = './data'
params.out_dir = './results'
params.rda_file = "${params.data_dir}/${params.study_id}.rda" // Assuming the RDA file is named after the study_id

// Define the process for loading the RDA file and extracting data
process LoadData {
    container 'nextflow-r-env:latest'

    input:
    path rda_file

    output:
    tuple val(params.study_id), path("expr_and_clin_data.Rds") into data_ch

    script:
    """
    Rscript load_data.R --rda_file ${params.rda_file} --output expr_and_clin_data.Rds
    """
}

// Define the process for running the gene survival analysis
process GeneSurvAnalysis {
    tag "$gene_idx"

    input:
    tuple val(study_id), path(data) from data_ch
    val gene_idx from 1..1000

    output:
    file("results_${gene_idx}.txt") into results_files

    script:
    """
    Rscript gene_surv_analysis.R --data ${data} --gene_idx $gene_idx --study_id $study_id \
                                  --out_file results_${gene_idx}.txt
    """
}

// Collate results into a single file
process CollateResults {
    input:
    file results from results_files.collect()

    output:
    file("combined_results.txt")

    """
    cat ${results} > combined_results.txt
    """
}

// Example R script for data loading (load_data.R)
"""
# load_data.R
args <- commandArgs(trailingOnly = TRUE)
rda_file <- args[1]
output_file <- args[2]

load(rda_file) # Load the RDA file
saveRDS(list(expr = expr, clin = clin), file = output_file) # Save expr and clin data to an RDS file
"""

// Example R script for gene analysis (gene_surv_analysis.R)
"""
# gene_surv_analysis.R
args <- commandArgs(trailingOnly = TRUE)
data_file <- args[1]
gene_idx <- as.integer(args[2])
study_id <- args[3]
out_file <- args[4]

data <- readRDS(data_file)
expr <- data$expr
clin <- data$clin

# Perform analysis...
# Save results...
"""

