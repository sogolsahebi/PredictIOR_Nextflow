#!/usr/bin/env Rscript
dir <- "./data/*.rda"

# Load a list of files of rda files
expr <- lapply(1:length(list.files(dir)), function(k){ 
load(file.path(dir, list.files(dir)[k])) 
dat_icb
})

# Assign study names to the list elements
study_icb <- substr(list.files(dir), 5, nchar(list.files(dir)) - 4)
names(expr) <- study_names

# Apply a function over the loaded datasets to perform survival analysis
cox_os <- lapply(1:length(expr), function(i) {
    study_parts <- unlist(strsplit(study_names[i], "__"))
    res <- geneSurvCont(
        dat.icb = expr[[i]],
        time.censor = 36,        // Censor time at 36 months
        missing.perc = 0.5,      // Exclude if more than 50% data is missing
        const.int = 0.001,       // A small constant to avoid division by zero or log(0)
        n.cutoff = 15,           // Minimum number of samples required
        feature = "CXCL9",       // Gene of interest
        study = study_names[i],  // Study ID
        surv.outcome = "OS",     // Analyze overall survival
        cancer.type = study_parts[2],  // Cancer type
        treatment = study_parts[3]     // Treatment type
    )

    # Remove rows with NA coefficients from the results
    res <- res[!is.na(res$Coef), ]
    res
})

# Combine the results for meta-analysis
meta_data <- do.call(rbind, cox_os)

# Perform meta-analysis using random-effects model
meta_analysis <- rma(yi = meta_data$HR, sei = meta_data$SE, method = "DL")

# Save the results to a CSV file
write.csv(meta_analysis, file = "ICB_Ravi__Lung__PD-L1_meta_analysis_os.csv", row.names = FALSE)
