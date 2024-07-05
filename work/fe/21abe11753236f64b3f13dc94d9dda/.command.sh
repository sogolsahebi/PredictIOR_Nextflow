#!/usr/bin/env Rscript
library(SummarizedExperiment)
library(PredictioR)
library(survival)
library(meta)

expr <- list()
for (dataset in dataset_names) {
    rda_path <- file.path("./data", paste0(dataset, ".rda"))
    load(rda_path)
    expr[[dataset]] <- get(dataset)
}

# Apply a function over the loaded datasets to perform survival analysis
assoc_res <- lapply(seq_along(expr), function(k) {
    geneSurvCont(
        dat.icb = expr[[k]],
        time.censor = 36,
        missing.perc = 0.5,
        const.int = 0.001,
        n.cutoff = 15,
        feature = "CXCL9",
        study = names(expr)[k],
        surv.outcome = 'OS',
        cancer.type = "c("Melanoma", "Pancreas", "Melanoma", "Bladder", "Melanoma", "Melanoma", "Kidney", "Melanoma")"[k],
        treatment = "c("PD-1/PD-L1", "PD-1/PD-L1", "PD-1/PD-L1", "PD-1/PD-L1", "CTLA4", "IO+combo", "PD-1/PD-L1", "CTLA4")"[k]
    )
})

# Combine results into one data frame and adjust for multiple testing
assoc_res <- do.call(rbind, assoc_res)
assoc_res$FDR <- p.adjust(assoc_res$Pval, method = "BH")

# Meta-analysis across datasets
res_meta_pancancer <- metafun(
    coef = assoc_res$Coef,
    se = assoc_res$SE,
    study = assoc_res$Study,
    pval = assoc_res$Pval,
    n = assoc_res$N,
    cancer.type = assoc_res$Cancer_type,
    treatment = assoc_res$Treatment,
    feature = "CXCL9",
    cancer.spec = FALSE,
    treatment.spec = FALSE
)

# Save the results to CSV files
write.csv(assoc_res, file = "Meta_analysis_os_CXCL9_association.csv", row.names = FALSE)
write.csv(res_meta_pancancer, file = "Meta_analysis_os_CXCL9_pancancer.csv", row.names = FALSE)
