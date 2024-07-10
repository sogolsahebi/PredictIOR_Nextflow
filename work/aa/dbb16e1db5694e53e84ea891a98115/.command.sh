#!/usr/bin/env Rscript
source('/R/load_libraries.R')

signature <- read.csv("signature_information.csv")
signature$Signature <- as.character(signature$signature)
signature$method <- as.character(signature$method)

dir_GeneSig <- "./SIG_data"
GeneSig_files <- list.files(path = dir_GeneSig, pattern = '*.rda', full.names = TRUE)
GeneSig_list <- lapply(GeneSig_files, function(file_path) {
    load(file_path)
    return(gsub(".rda", "", basename(file_path)))
})

load("ICB_small_Mariathasan.rda")

geneSig.score <- lapply(GeneSig_list, function(sig_name) {
    sig <- get(sig_name)

    method <- signature[signature$Signature == sig_name, "method"]
    if (method == "GSVA") {
        geneSig <- geneSigGSVA(...)
    } else if (method == "Weighted Mean") {
        geneSig <- geneSigMean(...)
    } else if (method == "ssGSEA") {
        geneSig <- geneSigssGSEA(...)
    } else if (method == "Specific Algorithm") {
        geneSig <- switch(sig_name,
                          "COX-IS_Bonavita" = geneSigCOX_IS(...),
                          "IPS_Charoentong" = geneSigIPS(...),
                          "PredictIO_Bareche" = geneSigPredictIO(...),
                          "IPRES_Hugo" = geneSigIPRES(...),
                          "PassON_Du" = geneSigPassON(...),
                          "IPSOV_Shen" = geneSigIPSOV(...))
    }

    if (length(geneSig) && !all(is.na(geneSig))) {
        return(geneSig)
    } else {
        return(rep(NA, length(sig)))
    }
})

geneSig.score <- do.call(rbind, geneSig.score)
write.csv(geneSig.score, file = "ICB_small_Mariathasan_GeneSigScore.csv", row.names = TRUE)
