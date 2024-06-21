#!/usr/bin/env Rscript
    library(SummarizedExperiment)
    library(PredictioR)
    library(survival)

    expr <- as.matrix(read.csv("ICB_small_Mariathasan_expr.csv", row.names = 1))
    clin <- read.csv("ICB_small_Mariathasan_clin.csv")

    # Print the first few rows of clinical data for debugging
    cat("Clinical data (first few rows):
")
    print(head(clin))

    # Ensure cancer_type is extracted correctly
    cancer_type <- names(table(clin$cancer_type)[table(clin$cancer_type) >= 15])
    cat("Cancer types with at least 15 samples:", cancer_type, "
")

    # Subset the clinical data to match the cancer type and sample IDs in the expression data
    relevant_samples <- clin$sample_id[clin$cancer_type %in% cancer_type]
    expr <- expr[, relevant_samples, drop=FALSE]
    clin <- clin[clin$sample_id %in% relevant_samples,]

    # Check the dimensions of the filtered expression matrix and clinical data
    cat("Filtered expression matrix dimensions:", dim(expr), "
")
    cat("Filtered clinical data dimensions:", dim(clin), "
")

    # Perform survival analysis
    status <- clin$event_occurred_os
    time <- clin$survival_time_os
    var <- as.numeric(scale(expr["CXCL9", ]))

    # Check lengths
    cat("Length of status:", length(status), "
")
    cat("Length of time:", length(time), "
")
    cat("Length of variable:", length(var), "
")

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
