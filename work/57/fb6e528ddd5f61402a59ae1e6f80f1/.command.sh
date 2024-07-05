#!/usr/bin/env Rscript
library(SummarizedExperiment)
library(PredictioR)
library(survival)
library(meta)

# Load all .rda files into the global environment
lapply(list.files(path = 'data', pattern = '*.rda', full.names = TRUE), function(file) {
    load(file, envir = .GlobalEnv)
})

# Create a list of the loaded objects using their actual names
loaded_objects <- ls(pattern = "^ICB_small_")
expr <- mget(loaded_objects, envir = .GlobalEnv)

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
                 cancer.type = "c("Melanoma", "Pancreas", "Melanoma", "Bladder", "Melanoma", "Melanoma", "Kidney", "Melanoma")"[k],
                 treatment = "c("PD-1/PD-L1", "PD-1/PD-L1", "PD-1/PD-L1", "PD-1/PD-L1", "CTLA4", "IO+combo", "PD-1/PD-L1", "CTLA4")"[k])
})

# Combine results into one data frame
assoc.res <- do.call(rbind, assoc.res)

# Adjust p-values for multiple comparisons using the Benjamini-Hochberg method
assoc.res$FDR <- p.adjust(assoc.res$Pval, method = "BH")

# meta-analysis for a gene across datasets
res_meta_pancancer <- metafun(coef = assoc.res$Coef, 
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
res_meta_pancancer <- data.frame(res_meta_pancancer)

# Save the results to a CSV file
write.csv(assoc.res, file = "Meta_analysis_os_CXCL9_assosiation.csv", row.names = FALSE)
write.csv(res_meta_pancancer, file = "Meta_analysis_os_CXCL9_pancancer.csv", row.names = FALSE)