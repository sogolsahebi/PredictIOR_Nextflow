#!/usr/bin/env Rscript
library(SummarizedExperiment)
library(PredictioR)
library(GSVA)

# Load ICB dataset
load(icb_data)
expr <- assay(ICB_small_Mariathasan)
clin <- as.data.frame(colData(ICB_small_Mariathasan))

# List to store signature scores
signature_scores_list <- list()

# Process each signature method
for (sig_name in params.signature_methods) {
    # Load signature file
    sig_file <- file.path(sig_data_dir, paste0(sig_name, ".rda"))
    load(sig_file)
    sig <- get("sig", envir = .GlobalEnv)

    # Determine method for processing
    method <- switch(sig_name,
                     ADO_Sidders = 'GSVA',
                     APM_Thompson = 'Weighted Mean',
                     APM_Wang = 'GSVA',
                     default = 'GSVA'  # Default method if not specified
    )

    # Process signature based on method
    if (method == 'GSVA') {
        geneSig <- geneSigGSVA(dat.icb = expr,
                               sig = sig,
                               sig.name = sig_name,
                               missing.perc = 0.5,
                               const.int = 0.001,
                               n.cutoff = 15,
                               sig.perc = 0.8,
                               study = "ICB_small_Mariathasan")
    } else if (method == 'Weighted Mean') {
        geneSig <- geneSigMean(dat.icb = expr,
                               sig = sig,
                               sig.name = sig_name,
                               missing.perc = 0.5,
                               const.int = 0.001,
                               n.cutoff = 15,
                               sig.perc = 0.8,
                               study = "ICB_small_Mariathasan")
    }
    # Add more methods as needed

    # Store signature scores
    signature_scores_list[[sig_name]] <- geneSig
}

# Save all signature scores as CSV
write.csv(do.call(rbind, signature_scores_list), file = "./output/ICB_small_Mariathasan_SignatureScores.csv", row.names = TRUE)
