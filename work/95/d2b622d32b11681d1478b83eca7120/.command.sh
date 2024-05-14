#!/usr/bin/env Rscript
library(SummarizedExperiment)

load("ICB_Ravi__Lung__PD-L1.rda")
expr <- assay(dat_icb)
clin <- as.data.frame(colData(dat_icb))

# Source the R script containing the geneSurvCont function and potentially others
source('/R_scripts/getGeneAssociation.R')
source('/R_scripts/getHR.R')

# Perform analysis
results <- lapply(1:100, function(i) {
  geneSurvCont(
    dat.icb = expr,  # Make sure 'expr' is defined in the loaded .rda or in getGeneAssociation.R
    clin = clin,     # Make sure 'clin' is defined in the loaded .rda or in getGeneAssociation.R
    time.censor = 36,
    missing.perc = 0.5,
    const.int = 0.001,
    n.cutoff = 15,
    feature = rownames(expr)[i],
    study = paste(study_id_parts[1], study_id_parts[2], study_id_parts[3], sep=''),
    surv.outcome = 'OS',
    cancer.type = study_id_parts[2],
    treatment = study_id_parts[3]
  )
})

# Save the results
write.csv(results, file = output_file_path, row.names = FALSE)
