#!/usr/bin/env Rscript
source('./load_libraries.R')

# Load all .rda files into the global environment
lapply(list.files(path = 'ICB_data', pattern = '*.rda', full.names = TRUE), function(file) {
    load(file, envir = .GlobalEnv)
})

# Create a list of the loaded objects using their actual names
loaded_objects <- ls(pattern = "^ICB_")
expr <- mget(loaded_objects, envir = .GlobalEnv)

# Define the cancer types and treatment types vectors
cancer_types <- fromJSON('["Melanoma", "Pancreas", "Melanoma", "Bladder", "Melanoma", "Melanoma", "Kidney", "Melanoma"]')
treatment_types <- fromJSON('["PD-1/PD-L1", "PD-1/PD-L1", "PD-1/PD-L1", "PD-1/PD-L1", "CTLA4", "IO+combo", "PD-1/PD-L1", "CTLA4"]')

# Apply a function over the loaded datasets to perform survival analysis
assoc.res <- lapply(1:length(expr), function(k){
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
})

# Combine results into one data frame
assoc.res <- do.call(rbind, assoc.res)

# Adjust p-values for multiple comparisons using the Benjamini-Hochberg method
assoc.res$FDR <- p.adjust(assoc.res$Pval, method = "BH")

# Treatment-specific meta-analysis for a gene across datasets
res_meta_percancer <- metaPerCanfun(
    coef = assoc.res$Coef, 
    se = assoc.res$SE,
    study  = assoc.res$Study, 
    pval = assoc.res$Pval, 
    n = assoc.res$N,
    cancer.type = assoc.res$Cancer_type,
    treatment = assoc.res$Treatment,
    feature = "CXCL9", 
    cancer.spec = TRUE
)

res_meta_percancer <- data.frame(rbind(res_meta_percancer$Melanoma$meta_summery, res_meta_percancer$Other$meta_summery))

# Save the results to a CSV file
write.csv(res_meta_percancer, file = "Meta_analysis_os_CXCL9_percancer.csv", row.names = FALSE)
