#!/usr/bin/env Rscript

# Callig the functions
source('/R_scripts/getGeneAssociation.R')
source('/R_scripts/getHR.R')

# Load each RDA file and extract data
list_rda <- lapply(list.files(path = 'data', pattern = '.rda', full.names = TRUE), function(file) {
    load(file)
    dat_icb
})

# Assign study names to the list elements
study_names <- substr(list.files('data'), 5, nchar(list.files('data')) - 4)
names(list_rda) <- study_names


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

# Combine the results for meta-analysis
meta_data <- do.call(rbind, cox_os)

# Perform meta-analysis using random-effects model
meta_analysis <- rma(yi = meta_data$HR, sei = meta_data$SE, method = "DL")

# Save the results to a CSV file
write.csv(meta_analysis, file = "ICB_Ravi__Lung__PD-L1_meta_analysis_os.csv", row.names = FALSE)
