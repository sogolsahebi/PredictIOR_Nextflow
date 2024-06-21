#!/usr/bin/env Rscript
library(SummarizedExperiment)
library(PredictioR)
library(survival)
library(meta)

list_rda <- lapply(list.files(path = 'data', pattern = '*.rda', full.names = TRUE), function(file) {
    load(file)
})

expr <- list('ICB_Liu' = ICB_small_Liu, 'ICB_Padron' = ICB_small_Padron, 'ICB_Hugo' = ICB_small_Hugo, 
             'ICB_Mariathasan' = ICB_small_Mariathasan, 'ICB_Nathanson' = ICB_small_Nathanson, 
             'ICB_Riaz' = ICB_small_Riaz, 'ICB_Miao' = ICB_small_Miao, 'ICB_Van_Allen' = ICB_small_Van_Allen)

cancer_type <- c('Melanoma', 'Pancreas', 'Melanoma', 'Bladder', 'Melanoma', 'Melanoma', 'Kidney', 'Melanoma')
treatment_type <- c('PD-1/PD-L1', 'PD-1/PD-L1', 'PD-1/PD-L1', 'PD-1/PD-L1', 'CTLA4', 'IO+combo', 'PD-1/PD-L1', 'CTLA4')

# Apply a function over the loaded datasets to perform survival analysis
assoc.res <- lapply(seq_along(expr), function(k){
    geneSurvCont(dat.icb = expr[[k]],
                 time.censor = 36,
                 missing.perc = 0.5,
                 const.int = 0.001,
                 n.cutoff = 15,
                 feature = "null",
                 study = names(expr)[k],
                 surv.outcome = 'OS',
                 cancer.type = cancer_type[k],
                 treatment = treatment_type[k])
})

# Combine results into one data frame
assoc.res <- do.call(rbind, assoc.res)

# Adjust p-values for multiple comparisons using the Benjamini-Hochberg method
assoc.res$FDR <- p.adjust(assoc.res$Pval, method = "BH")

res_meta_percancer <- metaPerCanfun(coef = assoc.res$Coef, 
                                    se = assoc.res$SE,
                                    study  = assoc.res$Study, 
                                    pval = assoc.res$Pval, 
                                    n = assoc.res$N,
                                    cancer.type = assoc.res$Cancer_type,
                                    treatment = assoc.res$Treatment,
                                    feature = "null", 
                                    cancer.spec = TRUE)

# Check if the `Other` category exists in the output structure
if (!is.null(res_meta_percancer$Other)) {
    res_meta_percancer  <- data.frame(rbind(res_meta_percancer$Melanoma$meta_summery, res_meta_percancer$Other$meta_summery))
} else {
    res_meta_percancer  <- data.frame(res_meta_percancer$Melanoma$meta_summery)
}

# Save the results to a CSV file
write.csv(res_meta_percancer, file = "ICB_small_Mariathasan_meta_analysis_os_null_percancer.csv", row.names = FALSE)
