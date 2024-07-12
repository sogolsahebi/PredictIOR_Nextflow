#!/usr/bin/env Rscript
source('/R/load_libraries.R')

signature <- read.csv("signature_information.csv")
signature$Signature <- as.character(signature$Signature)
signature$method <- as.character(signature$method)

dir_GeneSig <- './SIG_data'
GeneSig_list <- list.files(path = './SIG_data', pattern = '*.rda', full.names = TRUE)

load("ICB_small_Padron.rda")

computeGeneSignature <- function(method, sig_name, sig) {
    if (method == "GSVA") {
        geneSig <- geneSigGSVA(dat.icb = ICB_small_Padron, sig = sig, sig.name = sig_name, missing.perc = 0.5, const.int = 0.001, n.cutoff = 15, sig.perc = 0.8, study = "ICB_small_Padron")
        if (sum(!is.na(geneSig)) > 0) {
            geneSig <- geneSig[1,]
        }
    } else if (method == "Weighted Mean") {
        geneSig <- geneSigMean(dat.icb = ICB_small_Padron, sig = sig, sig.name = sig_name, missing.perc = 0.5, const.int = 0.001, n.cutoff = 15, sig.perc = 0.8, study = "ICB_small_Padron")
    } else if (method == "ssGSEA") {
        geneSig <- geneSigssGSEA(dat.icb = ICB_small_Padron, sig = sig, sig.name = sig_name, missing.perc = 0.5, const.int = 0.001, n.cutoff = 15, sig.perc = 0.8, study = "ICB_small_Padron")
        if (sum(!is.na(geneSig)) > 0) {
            geneSig <- geneSig[1,]
        }
    } else if (method == "Specific Algorithm") {
        geneSig <- geneSigPassON(dat.icb = ICB_small_Padron, sig = sig, sig.name = sig_name, missing.perc = 0.5, const.int = 0.001, n.cutoff = 15, sig.perc = 0.8, study = "ICB_small_Padron")
    }

    if (sum(!is.na(geneSig)) > 0) {
        return(geneSig)
    } else {
        return(rep(NA, ncol(ICB_small_Padron)))
    }
}

geneSig.score <- lapply(1:length(GeneSig_list), function(i) {
    load(GeneSig_list[i])
    sig_name <- substr(basename(GeneSig_list[i]), 1, nchar(basename(GeneSig_list[i])) - 4)
    method <- signature[signature$Signature == sig_name, "method"]
    computeGeneSignature(method, sig_name, sig)
})

geneSig.score <- do.call(rbind, geneSig.score)
rownames(geneSig.score) <- substr(basename(GeneSig_list), 1, nchar(basename(GeneSig_list)) - 4)
remove <- which(is.na(rowSums(geneSig.score)))
if (length(remove) > 0) {
    geneSig.score <- geneSig.score[-remove, ]
}
write.csv(geneSig.score, file = "ICB_small_Padron_GeneSigScore.csv", row.names = FALSE)
