#!/usr/bin/env Rscript

# Load necessary libraries
library(SummarizedExperiment)

# Command-line arguments are automatically passed by Nextflow
rda_file_path <- 'ICB_Ravi__Lung__PD-L1.rda'
output_file_path <- 'ICB_Ravi__Lung__PD-L1_results.csv'

# Load the data from the specified .rda file
load(rda_file_path)

# Extract study ID and its parts
study_id <- gsub('.rda', '', basename(rda_file_path))
study_id_parts <- unlist(strsplit(study_id, '__'))

# Source the R script containing the geneSurvCont function and potentially others
source('/R_scripts/getGeneAssociation.R')
source('/R_scripts/getHR.R')

load("ICB_Ravi__Lung__PD-L1.rda")

# Perform analysis
results <- lapply(1:100, function(i) {
  geneSurvCont(
    dat.icb = dat_icb,  # Make sure 'expr' is defined in the loaded .rda or in getGeneAssociation.R
    clin = NULL,     # Make sure 'clin' is defined in the loaded .rda or in getGeneAssociation.R
    time.censor = 36,
    missing.perc = 0.5,
    const.int = 0.001,
    n.cutoff = 15,
    feature = rownames(expr)[i],
    study = paste(study_id_parts[1], study_id_parts[2], study_id_parts[3], sep='__'),
    surv.outcome = 'OS',
    cancer.type = study_id_parts[2],
    treatment = study_id_parts[3]
  )
})

# Save the results
write.csv(results, file = output_file_path, row.names = FALSE)
