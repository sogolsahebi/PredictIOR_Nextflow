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
cancer_types <- fromJSON('["Melanoma", "Pancreas", "Melanoma", "Bladder",  "Melanoma", "Melanoma", "Kidney",   "Melanoma"]')
treatment_types <- fromJSON('["PD-1/PD-L1", "PD-1/PD-L1", "PD-1/PD-L1", "PD-1/PD-L1", "CTLA4", "IO+combo", "PD-1/PD-L1", "CTLA4"]')

# Apply a function over the loaded datasets to perform survival or response analysis
assoc.res <- lapply(1:length(expr), function(k){
    if ('OS' == 'OS' ) {
        geneSurvCont(
            dat.icb = expr[[k]],
            time.censor = 36,
            missing.perc = 0.5,
            const.int = 0.001,
            n.cutoff = 15,
            feature = "CXCL9",
            study = names(expr)[k],
            surv.outcome = 'OS',
            cancer.type = cancer_types[k],
            treatment = treatment_types[k]
        )
    } else if ('OS' == 'PFS') {
        geneSurvCont(
            dat.icb = expr[[k]],
            time.censor = 24,
            missing.perc = 0.5,
            const.int = 0.001,
            n.cutoff = 15,
            feature = "CXCL9",
            study = names(expr)[k],
            surv.outcome = 'OS',
            cancer.type = cancer_types[k],
            treatment = treatment_types[k]
        )
    } else if ('OS' == 'Response') {
        geneLogReg(
            dat.icb = expr[[k]],
            missing.perc = 0.5,
            const.int = 0.001,
            n.cutoff = 15,
            feature = "CXCL9",
            study = names(expr)[k],
            n0.cutoff = 3, 
            n1.cutoff = 3,  
            cancer.type = cancer_types[k],
            treatment = treatment_types[k]
        )
    }
})

# Combine results into one data frame
assoc.res <- do.call(rbind, assoc.res)

# Adjust p-values for multiple comparisons using the Benjamini-Hochberg method
assoc.res$FDR <- p.adjust(assoc.res$Pval, method = "BH")

# Meta-analysis for a gene across datasets

# Treatment-specific meta-analysis for a gene across datasets
 res_meta_pancancer <- metafun(
    coef = assoc.res$Coef, 
    se = assoc.res$SE,
    study  = assoc.res$Study, 
    pval = assoc.res$Pval, 
    n = assoc.res$N,
    cancer.type = assoc.res$Cancer_type,
    treatment = assoc.res$Treatment,
    feature = "CXCL9", 
    cancer.spec = FALSE, 
    treatment.spec = FALSE
)

# Save the results to a CSV file
write.csv(data.frame(res_meta_pancancer), file = "Meta_analysis_OS_CXCL9_pancancer.csv", row.names = FALSE)
