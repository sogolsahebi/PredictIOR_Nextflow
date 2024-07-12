#!/usr/bin/env Rscript
source('/R/load_libraries.R')

# Load all .rda files 
lapply(list.files(path = 'ICB_data', pattern = '*.rda', full.names = TRUE), function(file) { 
    load(file, envir = .GlobalEnv)
})

# Create a list of the loaded objects using their actual names
loaded_objects <- basename(list.files(path = 'ICB_data', pattern = '*.rda', full.names = TRUE))
loaded_objects <- substr(loaded_objects, 1, nchar(loaded_objects) - 4) # remove .rda from name
expr <- mget(loaded_objects, envir = .GlobalEnv)

# Define the cancer types and treatment types vectors
cancer_types <- fromJSON('["Melanoma", "Pancreas", "Melanoma", "Bladder", "Melanoma", "Melanoma", "Kidney", "Melanoma"]')
treatment_types <- fromJSON('["PD-1/PD-L1", "PD-1/PD-L1", "PD-1/PD-L1", "PD-1/PD-L1", "CTLA4", "IO+combo", "PD-1/PD-L1", "CTLA4"]')

# Apply a function over the loaded datasets to perform survival or response analysis
assoc.res <- lapply(1:length(expr), function(k){
    print(paste("Processing dataset:", names(expr)[k])) # Debug
    if ('Response' == 'OS') {
        return(geneSurvCont(dat.icb = expr[[k]], time.censor = 36, missing.perc = 0.5, const.int = 0.001, n.cutoff = 15, feature = "CXCL9", study = names(expr)[k], surv.outcome = 'Response', cancer.type = cancer_types[k], treatment = treatment_types[k]))
    } else if ('Response' == 'PFS') {
        return(geneSurvCont(dat.icb = expr[[k]], time.censor = 24, missing.perc = 0.5, const.int = 0.001, n.cutoff = 15, feature = "CXCL9", study = names(expr)[k], surv.outcome = 'Response', cancer.type = cancer_types[k], treatment = treatment_types[k]))
    } else if ('Response' == 'Response') {
        return(geneLogReg(dat.icb = expr[[k]], missing.perc = 0.5, const.int = 0.001, n.cutoff = 15, feature = "CXCL9", study = names(expr)[k], n0.cutoff = 3, n1.cutoff = 3, cancer.type = cancer_types[k], treatment = treatment_types[k]))
    } else {
        print(paste("You can only use 'OS', 'PFS', or 'Response' as 'Response' input."))
        return(NULL)
    }
})

# Combine results into one data frame
assoc.res <- do.call(rbind, assoc.res)

# Adjust p-values for multiple comparisons using the Benjamini-Hochberg method
assoc.res$FDR <- p.adjust(assoc.res$Pval, method = "BH")

# Treatment-specific meta-analysis for a gene across datasets
res_meta_percancer <- metaPerCanfun(coef = assoc.res$Coef, se = assoc.res$SE, study = assoc.res$Study, pval = assoc.res$Pval, n = assoc.res$N, cancer.type = assoc.res$Cancer_type, treatment = assoc.res$Treatment, feature = "CXCL9", cancer.spec = TRUE)

# Combine all meta_summery results into a single data frame
meta_summery_combined <- do.call(rbind, lapply(res_meta_percancer, function(x) x$meta_summery))
write.csv(meta_summery_combined, file = "Meta_analysis_Response_CXCL9_percancer.csv", row.names = FALSE)