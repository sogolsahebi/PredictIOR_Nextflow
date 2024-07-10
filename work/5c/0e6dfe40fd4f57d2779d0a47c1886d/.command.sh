#!/usr/bin/env Rscript
library(PredictioR)
library(survival)
library(GSVA)
library(dplyr)
library(data.table)
library(MultiAssayExperiment)
library(Hmisc)

# Load signature information
signature <- read.csv("signature_information.csv")
signature$Signature <- as.character(signature$signature)
signature$method <- as.character(signature$method)

# Load gene signature files
dir_GeneSig <- './SIG_data'
GeneSig_list <- list.files(path = dir_GeneSig, pattern = '*.rda', full.names = TRUE)

# Load ICB data
load("ICB_small_Mariathasan.rda")

# Function to compute signature score based on method
compute_signature_score <- function(method, dat_icb, sig, sig_name) {
    if (method == "GSVA") {
        geneSig <- geneSigGSVA(dat.icb = dat_icb, sig = sig, sig.name = sig_name, missing.perc = 0.5, const.int = 0.001, n.cutoff = 15, sig.perc = 0.8, study = "ICB_small_Mariathasan")
    } else if (method == "Weighted Mean") {
        geneSig <- geneSigMean(dat.icb = dat_icb, sig = sig, sig.name = sig_name, missing.perc = 0.5, const.int = 0.001, n.cutoff = 15, sig.perc = 0.8, study = "ICB_small_Mariathasan")
    } else if (method == "ssGSEA") {
        geneSig <- geneSigssGSEA(dat.icb = dat_icb, sig = sig, sig.name = sig_name, missing.perc = 0.5, const.int = 0.001, n.cutoff = 15, sig.perc = 0.8, study = "ICB_small_Mariathasan")
    } else if (method == "Specific Algorithm") {
        switch(sig_name,
            "COX-IS_Bonavita" = geneSig <- geneSigCOX_IS(dat.icb = dat_icb, sig = sig, sig.name = sig_name, missing.perc = 0.5, const.int = 0.001, n.cutoff = 15, sig.perc = 0.8, study = "ICB_small_Mariathasan"),
            "IPS_Charoentong" = geneSig <- geneSigIPS(dat.icb = dat_icb, sig = sig, sig.name = sig_name, missing.perc = 0.5, const.int = 0.001, n.cutoff = 15, sig.perc = 0.8, study = "ICB_small_Mariathasan"),
            "PredictIO_Bareche" = geneSig <- geneSigPredictIO(dat.icb = dat_icb, sig = sig, sig.name = sig_name, missing.perc = 0.5, const.int = 0.001, n.cutoff = 15, sig.perc = 0.8, study = "ICB_small_Mariathasan"),
            "IPRES_Hugo" = geneSig <- geneSigIPRES(dat.icb = dat_icb, sig = sig, sig.name = sig_name, missing.perc = 0.5, const.int = 0.001, n.cutoff = 15, sig.perc = 0.8, study = "ICB_small_Mariathasan"),
            "PassON_Du" = geneSig <- geneSigPassON(dat.icb = dat_icb, sig = sig, sig.name = sig_name, missing.perc = 0.5, const.int = 0.001, n.cutoff = 15, sig.perc = 0.8, study = "ICB_small_Mariathasan"),
            "IPSOV_Shen" = geneSig <- geneSigIPSOV(dat.icb = dat_icb, sig = sig, sig.name = sig_name, missing.perc = 0.5, const.int = 0.001, n.cutoff = 15, sig.perc = 0.8, study = "ICB_small_Mariathasan"),
            stop("Unsupported specific algorithm: ", sig_name)
        )
    } else {
        stop("Unsupported method: ", method)
    }
    if (sum(!is.na(geneSig)) == 0) {
        geneSig <- rep(NA, ncol(dat_icb))
    }
    geneSig
}

geneSig.score <- lapply(GeneSig_list, function(gene_sig_file) {
    load(gene_sig_file)
    sig_name <- substr(basename(gene_sig_file), 1, nchar(basename(gene_sig_file)) - 4)
    method <- signature[signature$Signature == sig_name, "method"]
    geneSig <- compute_signature_score(method, dat_icb = expr, sig = sig, sig_name = sig_name)
    geneSig
})

geneSig.score <- do.call(rbind, geneSig.score)
rownames(geneSig.score) <- substr(basename(GeneSig_list), 1, nchar(basename(GeneSig_list)) - 4)

remove <- which(is.na(rowSums(geneSig.score)))
if (length(remove) > 0) {
    geneSig.score <- geneSig.score[-remove, ]
}
write.csv(geneSig.score, file = "ICB_small_Mariathasan_GeneSigScore.csv", row.names = FALSE)
