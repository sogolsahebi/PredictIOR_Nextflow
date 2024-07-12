#!/usr/bin/env Rscript
source('/R/load_libraries.R')

signature <- read.csv("signature_information.csv")
signature$Signature <- as.character(signature$signature)
signature$method <- as.character(signature$method)

dir_GeneSig <- './SIG_data'
GeneSig_list <- list.files(path = './SIG_data', pattern = '*.rda', full.names = TRUE)

load("ICB_small_Miao.rda")

geneSig.score <- lapply(1:length(GeneSig_list), function(i) {
    load(GeneSig_list[i])
    sig_name <- substr(basename(GeneSig_list[i]), 1, nchar(basename(GeneSig_list[i])) - 4)


    method <- signature[signature$Signature == sig_name, "method"]

    if (signature[signature$Signature == sig_name, "method"] == "GSVA") {
        geneSig <- geneSigGSVA(dat.icb = ICB_small_Miao, sig = sig, sig.name = sig_name, missing.perc = 0.5, const.int = 0.001, n.cutoff = 15, sig.perc = 0.8, study = "ICB_small_Miao")
        if (sum(!is.na(geneSig)) > 0) {
            geneSig <- geneSig[1,]
        }
    } else if (signature[signature$Signature == sig_name, "method"] == "Weighted Mean") {
        geneSig <- geneSigMean(dat.icb = ICB_small_Miao, sig = sig, sig.name = sig_name, missing.perc = 0.5, const.int = 0.001, n.cutoff = 15, sig.perc = 0.8, study = "ICB_small_Miao")
    } else if (signature[signature$Signature == sig_name, "method"] == "ssGSEA") {
        geneSig <- geneSigssGSEA(dat.icb = ICB_small_Miao, sig = sig, sig.name = sig_name, missing.perc = 0.5, const.int = 0.001, n.cutoff = 15, sig.perc = 0.8, study = "ICB_small_Miao")
        if (sum(!is.na(geneSig)) > 0) {
            geneSig <- geneSig[1,]
        }
    } else if (signature[signature$Signature == sig_name, "method"] == "Specific Algorithm" & sig_name == "COX-IS_Bonavita") {
        geneSig <- geneSigCOX_IS(dat.icb = ICB_small_Miao, sig = sig, sig.name = signature$Signature[i], missing.perc = 0.5, const.int = 0.001, n.cutoff = 15, sig.perc = 0.8, study = "ICB_small_Miao")
    } else if (signature[signature$Signature == sig_name, "method"] == "Specific Algorithm" & sig_name == "IPS_Charoentong") {
        geneSig <- geneSigIPS(dat.icb = ICB_small_Miao, sig = sig, sig.name = signature$Signature[i], missing.perc = 0.5, const.int = 0.001, n.cutoff = 15, sig.perc = 0.8, study = "ICB_small_Miao")
    } else if (signature[signature$Signature == sig_name, "method"] == "Specific Algorithm" & sig_name == "PredictIO_Bareche") {
        geneSig <- geneSigPredictIO(dat.icb = ICB_small_Miao, sig = sig, sig.name = signature$Signature[i], missing.perc = 0.5, const.int = 0.001, n.cutoff = 15, sig.perc = 0.8, study = "ICB_small_Miao")
    } else if (signature[signature$Signature == sig_name, "method"] == "Specific Algorithm" & sig_name == "IPRES_Hugo") {
        geneSig <- geneSigIPRES(dat.icb = ICB_small_Miao, sig = sig, sig.name = signature$Signature[i], missing.perc = 0.5, const.int = 0.001, n.cutoff = 15, sig.perc = 0.8, study = "ICB_small_Miao")
    } else if (signature[signature$Signature == sig_name, "method"] == "Specific Algorithm" & sig_name == "PassON_Du") {
        geneSig <- geneSigPassON(dat.icb = ICB_small_Miao, sig = sig, sig.name = signature$Signature[i], missing.perc = 0.5, const.int = 0.001, n.cutoff = 15, sig.perc = 0.8, study = "ICB_small_Miao")
    } else if (signature[signature$Signature == sig_name, "method"] == "Specific Algorithm" & sig_name == "IPSOV_Shen") {
        geneSig <- geneSigIPSOV(dat.icb = ICB_small_Miao, sig = sig, sig.name = signature$Signature[i], missing.perc = 0.5, const.int = 0.001, n.cutoff = 15, sig.perc = 0.8, study = "ICB_small_Miao")
    }

    if (sum(!is.na(geneSig)) > 0) {
        geneSig <- geneSig
    } else {
        geneSig <- rep(NA, ncol(ICB_small_Miao))
    }

    geneSig
})

geneSig.score <- do.call(rbind, geneSig.score)
rownames(geneSig.score) <- substr(basename(GeneSig_list), 1, nchar(basename(GeneSig_list)) - 4)

remove <- which(is.na(rowSums(geneSig.score)))
if (length(remove) > 0) {
    geneSig.score <- geneSig.score[-remove, ]
}
write.csv(geneSig.score, file = "ICB_small_Miao_GeneSigScore.csv", row.names = TRUE)
