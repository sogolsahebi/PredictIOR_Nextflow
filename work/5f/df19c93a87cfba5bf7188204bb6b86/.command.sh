#!/usr/bin/env Rscript
source('/R/load_libraries.R')

# Load all cox results
cox_results <- lapply(ICB_small_Mariathasan_Response.csv, function(file) {
    read.csv(file)
})

# Combine results into one data frame
assoc.res <- do.call(rbind, cox_results)

# Meta-analysis for a gene across datasets
res_meta <- metaPerTreatmentfun(
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
res_meta <- data.frame(rbind(res_meta$`PD-1/PD-L1`$meta_summery, res_meta$Other$meta_summery))

write.csv(data.frame(res_meta), file = "Meta_analysis_Response_CXCL9_pertreatment.csv", row.names = FALSE)
