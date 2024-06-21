#!/usr/bin/env Rscript
    library(SummarizedExperiment)
    library(PredictioR)

    expr <- as.matrix(read.csv("ICB_small_Mariathasan_expr.csv", row.names = 1))
    clin <- read.csv("ICB_small_Mariathasan_clin.csv")

    # Ensure cancer_type is extracted correctly
    cancer_type <- names(table(clin$cancer_type)[table(clin$cancer_type) >= 15])
    
    # Perform survival analysis
    status <- clin$event_occurred_os[clin$cancer_type %in% cancer_type]
    time <- clin$survival_time_os[clin$cancer_type %in% cancer_type]
    var <- as.numeric(scale(expr["CXCL9", clin$cancer_type %in% cancer_type]))
    
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
