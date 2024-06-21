#!/usr/bin/env Rscript
library(PredictioR)

# Load each RDA file and extract data
list_rda <- lapply(list.files(path = 'data', pattern = '*.rda', full.names = TRUE), function(file) {
    load(file)
})

cancer_type <- c('Melanoma', 'Pancreas', 'Melanoma', 'Bladder', 'Melanoma', 'Melanoma', 'Kidney', 'Melanoma')
treatment_type <- c('PD-1/PD-L1', 'PD-1/PD-L1', 'PD-1/PD-L1', 'PD-1/PD-L1', 'CTLA4', 'IO+combo', 'PD-1/PD-L1', 'CTLA4')

# Apply a function over the loaded datasets to perform survival analysis
assoc.res <- lapply(1:length(list_rda), function(k){
    geneSurvCont(dat.icb = list_rda[[k]],
                 time.censor = 36,
                 missing.perc = 0.5,
                 const.int = 0.001,
                 n.cutoff = 15,
                 feature = CXCL9,
                 study = names(list_rda)[k],
                 surv.outcome = 'OS',
                 cancer.type = cancer_type[k],
                 treatment = treatment_type[k])
})

# Combine results into one data frame
assoc.res <- do.call(rbind, assoc.res)

# Adjust p-values for multiple comparisons using the Benjamini-Hochberg method
assoc.res$FDR <- p.adjust(assoc.res$Pval, method = "BH")

# Save the results to a CSV file
write.csv(assoc.res, file = "ICB_small_Mariathasan_meta_analysis_os.csv", row.names = FALSE)
