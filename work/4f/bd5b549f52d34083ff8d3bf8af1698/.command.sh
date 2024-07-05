#!/usr/bin/env Rscript
library(SummarizedExperiment)
library(PredictioR)
library(survival)

list_rda <- lapply(list.files(path = 'data', pattern = '*.rda', full.names = TRUE), function(file) {
    load(file)
})

expr <- list('ICB_Liu' = ICB_small_Liu, 'ICB_Padron' = ICB_small_Padron, 'ICB_Hugo' = ICB_small_Hugo, 
          'ICB_Mariathasan' = ICB_small_Mariathasan, 'ICB_Nathanson' = ICB_small_Nathanson, 
          'ICB_Riaz' = ICB_small_Riaz, 'ICB_Miao' = ICB_small_Miao, 'ICB_Van_Allen' = ICB_small_Van_Allen)

cancer_type <- c('Melanoma', 'Pancreas', 'Melanoma', 'Bladder', 'Melanoma', 'Melanoma', 'Kidney', 'Melanoma')
treatment_type <- c('PD-1/PD-L1', 'PD-1/PD-L1', 'PD-1/PD-L1', 'PD-1/PD-L1', 'CTLA4', 'IO+combo', 'PD-1/PD-L1', 'CTLA4')

# Apply a function over the loaded datasets to perform survival analysis
assoc.res <- lapply(1:length(expr), function(k){
    geneSurvCont(dat.icb = expr[[k]],
                 time.censor = 36,
                 missing.perc = 0.5,
                 const.int = 0.001,
                 n.cutoff = 15,
                 feature = "CXCL9",
                 study = names(expr)[k],
                 surv.outcome = 'OS',
                 cancer.type = cancer_type[k],
                 treatment = treatment_type[k])
})

# Combine results into one data frame
assoc.res <- do.call(rbind, assoc.res)

# Adjust p-values for multiple comparisons using the Benjamini-Hochberg method
assoc.res$FDR <- p.adjust(assoc.res$Pval, method = "BH")


# meta-analysis for a gene across datasets
res_meta <- metafun(coef = assoc.res$Coef, 
                    se = assoc.res$SE,
                    study  = assoc.res$Study, 
                    pval = assoc.res$Pval, 
                    n = assoc.res$N,
                    cancer.type = assoc.res$Cancer_type,
                    treatment = assoc.res$Treatment,
                    feature = "CXCL9", 
                    cancer.spec = FALSE, 
                    treatment.spec = FALSE)

# meta-analysis results
res_meta <- data.frame(res_meta)

# Save the results to a CSV file
write.csv(res_meta, file = "ICB_small_Mariathasan_meta_analysis_os__CXCL9_across_datasets.csv", row.names = FALSE)