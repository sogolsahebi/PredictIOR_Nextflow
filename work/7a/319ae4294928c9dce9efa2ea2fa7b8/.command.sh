#!/usr/bin/env Rscript
source('/R/load_libraries.R')

# Load all .rda files 
lapply(list.files(path = 'ICB_data', pattern = '*.rda', full.names = TRUE), function(file) { load(file, envir = .GlobalEnv) })

# Create a list of the loaded objects using their actual names
loaded_objects <- basename(list.files(path = 'ICB_data', pattern = '*.rda', full.names = TRUE))
loaded_objects <- substr(loaded_objects, 1, nchar(loaded_objects) - 4) # remove .rda from name
expr <- mget(loaded_objects, envir = .GlobalEnv)

# Define the cancer types and treatment types vectors
cancer_types <- fromJSON('["Melanoma", "Pancreas", "Melanoma", "Bladder", "Melanoma", "Melanoma", "Kidney", "Melanoma"]')
treatment_types <- fromJSON('["PD-1/PD-L1", "PD-1/PD-L1", "PD-1/PD-L1", "PD-1/PD-L1", "CTLA4", "IO+combo", "PD-1/PD-L1", "CTLA4"]')


geneSig.score <- read.csv("signature_information.csv")

# Apply a function over the loaded datasets to perform survival or response analysis
assoc.res <- lapply(1:length(expr), function(k){
    if ('OS' == 'OS') {
    res.all <- lapply(1:nrow(geneSig.score), function(k) {
    sig_name <- rownames(geneSig.score)[k]
    geneSig_vector <- as.numeric(geneSig.score[k, ])
    geneSig_vector <- geneSig_vector[!is.na(geneSig_vector)]

    res <- geneSigSurvCont(
        dat.icb = expr[[k]],
        geneSig = geneSig_vector,
        time.censor = 36,
        n.cutoff = 15,
        study = names(expr)[k],
        surv.outcome = "OS",
        sig.name = sig_name,
        cancer.type = cancer_types[k],
        treatment = treatment_types[k]
    )

    res
})

} else if ('OS' == 'PFS') {

res <- geneSigSurvCont(
    dat.icb = expr[[k]],
    geneSig = geneSig_vector,  
    time.censor = 24,
    n.cutoff = 15,
    study = names(expr)[k],
    surv.outcome = "PFS",
    sig.name = sig_name,
    cancer.type = cancer_types[k],
    treatment = treatment_types[k]
)

res
})

} else if ('OS' == 'Response') {

res <- geneSigLogReg(dat.icb = expr[[k]],
                    geneSig = geneSig_vector,
                    n.cutoff = 10,
                    study =  names(expr)[k],
                    sig.name = sig_name,
                    n0.cutoff = 3, 
                    n1.cutoff = 3,
                    cancer.type = cancer_types[k],
                    treatment = treatment_types[k])

res
})

} else {
        print(paste("You can only use 'OS', 'PFS', or 'Response' as 'OS' input."))
        return(NULL)
}
})

assoc.res <- do.call(rbind, assoc.res)
assoc.res$FDR <- p.adjust(assoc.res$Pval, method="BH")
assoc.res <- assoc.res[order(assoc.res$FDR), ]

# Meta-analysis for a gene across datasets
res_meta_pancancer <- metafun(coef = assoc.assoc.res$Coef, se = assoc.res$SE, study = assoc.res$Study, pval = assoc.all$Pval, n = assoc.all$N, cancer.type = assoc.all$Cancer_type, treatment = assoc.all$Treatment, feature = "CXCL9", cancer.spec = FALSE, treatment.spec = FALSE)

# Save the results to a CSV file
write.csv(data.frame(res_meta_pancancer), file = "Meta_analysis_OS_CXCL9_pancancer.csv", row.names = FALSE)
