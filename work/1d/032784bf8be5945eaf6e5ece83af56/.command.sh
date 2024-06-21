#!/usr/bin/env Rscript
    library(SummarizedExperiment)
    library(PredictioR)
    library(survival)

    expr <- as.matrix(read.csv("ICB_small_Mariathasan_expr.csv", row.names = 1))
    clin <- read.csv("ICB_small_Mariathasan_clin.csv", row.names = 1)

    # Ensure cancer_type is extracted correctly
    cancer_type <- names(table(clin$cancer_type)[table(clin$cancer_type) >= 15])
    cat("Cancer types with at least 15 samples:", cancer_type, "\n")

    # Subset the clinical data to match the cancer type and sample IDs in the expression data
    relevant_samples <- rownames(clin)[clin$cancer_type %in% cancer_type]
    relevant_samples <- relevant_samples[relevant_samples %in% colnames(expr)]
    expr <- expr[, relevant_samples, drop=FALSE]
    clin <- clin[relevant_samples, , drop=FALSE]

    # Check the dimensions of the filtered expression matrix and clinical data
    cat("Filtered expression matrix dimensions:", dim(expr), "\n")
    cat("Filtered clinical data dimensions:", dim(clin), "\n")
event_occurred_os
    # Perform survival analysis
    status <- clin$event_occurred_os
    time <- clin$survival_time_os
    var <- as.numeric(scale(expr["CXCL9", ]))

    # Check lengths
    cat("Length of status:", length(status), "\n")
    cat("Length of time:", length(time), "\n")
    cat("Length of variable:", length(var), "\n")

    if(length(status) == length(time) && length(time) == length(var)) {
        cox_os <- survCont(
            status = status,
            time = time,
            time.censor = 36,
            var = var
        )

        # Convert to dataframe
        cox_os <- as.data.frame(t(cox_os))

        write.csv(cox_os, file = "ICB_small_Mariathasan_cox_os_CXCL9.csv", row.names = FALSE)
    } else {
        stop("Lengths of status, time, and variable do not match.")
    }
