#!/usr/bin/env Rscript
source('/R/load_libraries.R')

# Load all .rda files into the global environment
lapply(list.files(path = 'ICB_small_Hugo.rda ICB_small_Liu.rda ICB_small_Mariathasan.rda ICB_small_Miao.rda ICB_small_Nathanson.rda ICB_small_Padron.rda ICB_small_Riaz.rda ICB_small_Van_Allen.rda', pattern = '*.rda', full.names = TRUE), function(file) {
    load(file, envir = .GlobalEnv)
})

# Create a list of the loaded objects using their actual names
loaded_objects <- ls(pattern = "^ICB_")
print("Loaded objects:")
print(loaded_objects)

expr <- mget(loaded_objects, envir = .GlobalEnv)
print("Expression data:")
print(expr)

# Define the cancer types and treatment types vectors
cancer_types <- fromJSON('["Melanoma", "Pancreas", "Melanoma", "Bladder", "Melanoma", "Melanoma", "Kidney", "Melanoma"]')
treatment_types <- fromJSON('["PD-1/PD-L1", "PD-1/PD-L1", "PD-1/PD-L1", "PD-1/PD-L1", "CTLA4", "IO+combo", "PD-1/PD-L1", "CTLA4"]')

# Check the lengths of cancer_types and treatment_types
print(paste("Number of loaded objects:", length(expr)))
print(paste("Number of cancer types:", length(cancer_types)))
print(paste("Number of treatment types:", length(treatment_types)))

# Apply a function over the loaded datasets to perform survival or response analysis
assoc.res <- lapply(1:length(expr), function(k){
    print(paste("Processing:", k))
    print(paste("Cancer type:", cancer_types[k]))
    print(paste("Treatment type:", treatment_types[k]))
    if ('OS' == 'OS' || 'OS' == 'PFS') {
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
    } else if ('OS' == 'Response') {
        geneLogReg(
            dat.icb = expr[[k]],
            missing.perc = 0.5,
            const.int = 0.001,
            n.cutoff = 15,
            feature = "CXCL9",
            study = names(expr)[k],
            n0.cutoff = 5, 
            n1.cutoff = 5,
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
if ('OS' == 'OS' || 'OS' == 'PFS') {
    res_meta_pertreatment <- metaPerTreatfun(
        coef = assoc.res$Coef, 
        se = assoc.res$SE,
        study  = assoc.res$Study, 
        pval = assoc.res$Pval, 
        n = assoc.res$N,
        cancer.type = assoc.res$Cancer_type,
        treatment = assoc.res$Treatment,
        feature = "CXCL9", 
        treatment.spec = TRUE
    )
} else if ('OS' == 'Response') {
    res_meta_pertreatment <- metaLogRegfun(
        coef = assoc.res$Coef, 
        se = assoc.res$SE,
        study  = assoc.res$Study, 
        pval = assoc.res$Pval, 
        n = assoc.res$N,
        cancer.type = assoc.res$Cancer_type,
        treatment = assoc.res$Treatment,
        feature = "CXCL9", 
        treatment.spec = TRUE
    )
}

res_meta_pertreatment <- data.frame(rbind(res_meta_pertreatment$PD1$meta_summery, res_meta_pertreatment$CTLA4$meta_summery))

# Save the results to a CSV file
write.csv(res_meta_pertreatment, file = "Meta_analysis_OS_CXCL9_pertreatment.csv", row.names = FALSE)
