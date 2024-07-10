#!/usr/bin/env Rscript
source('/R/load_libraries.R')  // Ensure this script properly sets up your R environment

// Read signature information and prepare data
signature <- read.csv("signature_information.csv")
signature$Signature <- as.character(signature$signature)
signature$method <- as.character(signature$method)

// Load signature files from the specified directory
dir_GeneSig <- "./SIG_data"
GeneSig_files <- list.files(path = dir_GeneSig, pattern = '*.rda', full.names = TRUE)
GeneSig_list <- lapply(GeneSig_files, function(file_path) {
    load(file_path)
    return(gsub(".rda", "", basename(file_path)))
})

// Load immunotherapy data
load("ICB_small_Mariathasan.rda")

// Compute signature scores based on different methods
geneSig.score <- lapply(GeneSig_list, function(sig_name) {
    sig <- get(sig_name)
    method <- signature[signature$Signature == sig_name, "method"]

    // Compute the signature score based on the method specified
    if (method == "GSVA") {
        geneSig <- geneSigGSVA(dat.icb, sig, 0.2, 0.001, 15, 0.8, "ICB_small_Mariathasan")
    } else if (method == "Weighted Mean") {
        geneSig <- geneSigMean(dat.icb, sig, 0.2, 0.001, 15, 0.8, "ICB_small_Mariathasan")
    } else if (method == "ssGSEA") {
        geneSig <- geneSigssGSEA(dat.icb, sig, 0.2, 0.001, 15, 0.8, "ICB_small_Mariathasan")
    } else if (method == "Specific Algorithm") {
        // Handle specific algorithms based on the signature name
        geneSig <- switch(sig_name,
                          "COX-IS_Bonavita" = geneSigCOX_IS(dat.icb, sig, 0.2, 0.001, 15, 0.8, "ICB_small_Mariathasan"),
                          "IPS_Charoentong" = geneSigIPS(dat.icb, sig, 0.2, 0.001, 15, 0.8, "ICB_small_Mariathasan"),
                          "PredictIO_Bareche" = geneSigPredictIO(dat.icb, sig, 0.2, 0.001, 15, 0.8, "ICB_small_Mariathasan"),
                          "IPRES_Hugo" = geneSigIPRES(dat.icb, sig, 0.2, 0.001, 15, 0.8, "ICB_small_Mariathasan"),
                          "PassON_Du" = geneSigPassON(dat.icb, sig, 0.2, 0.001, 15, 0.8, "ICB_small_Mariathasan"),
                          "IPSOV_Shen" = geneSigIPSOV(dat.icb, sig, 0.2, 0.001, 15, 0.8, "ICB_small_Mariathasan")
                          )
    }

    // Return the signature score, handle NA values
    if (length(geneSig) && !all(is.na(geneSig))) {
        return(geneSig)
    } else {
        return(rep(NA, length(sig)))
    }
})

// Combine all scores into a single file
geneSig.score <- do.call(rbind, geneSig.score)
write.csv(geneSig.score, file = "ICB_small_Mariathasan_GeneSigScore.csv", row.names = TRUE)
