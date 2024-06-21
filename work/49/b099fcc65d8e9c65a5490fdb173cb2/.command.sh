#!/usr/bin/env Rscript

# Load prepared data
prepared_data <- readRDS('prepared_data.rds')
list_rda <- prepared_data[[1]]
study_names <- prepared_data[[2]]

# Apply a function over the loaded datasets to perform survival analysis
cox_os <- lapply(1:length(list_rda), function(i) {
    study_parts <- unlist(strsplit(study_names[i], "__"))
    res <- geneSurvCont(
        dat.icb = list_rda[[i]],
        time.censor = 36,
        missing.perc = 0.5,
        const.int = 0.001,
        n.cutoff = 15,
        feature = "CXCL9",
        study = study_names[i],
        surv.outcome = "OS",
        cancer.type = study_parts[2],
        treatment = study_parts[3]
    )

    # Remove rows with NA coefficients from the results
    res <- res[!is.na(res$Coef), ]
    res
})

# Combine results into one data frame
cox_os <- do.call(rbind, cox_os)

# Remove rows without a gene name
cox_os <- cox_os[!is.na(cox_os$Gene), ]

# Adjust p-values for multiple comparisons using the Benjamini-Hochberg method
cox_os$FDR <- p.adjust(cox_os$Pval, method = "BH")

# Save the results to a CSV file
write.csv(cox_os, file = "ICB_Ravi__Lung__PD-L1_meta_analysis_os.csv", row.names = FALSE)
