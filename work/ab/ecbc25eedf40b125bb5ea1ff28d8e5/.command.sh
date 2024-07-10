#!/usr/bin/env Rscript
source('/R/load_libraries.R')
signature <- read.csv("signature_information.csv")
signature$Signature <- as.character(signature$signature)
signature$method <- as.character(signature$method)

dir_GeneSig <- './SIG_data'
GeneSig_list <- list.files(path = './SIG_data', pattern = '*.rda', full.names = TRUE)

load("ICB_small_Mariathasan.rda")

geneSig.score <- lapply(1:length(GeneSig_list), function(i) {
    load(GeneSig_list[i])
    sig_name <- substr(basename(GeneSig_list[i]), 1, nchar(basename(GeneSig_list[i])) - 4)
    method <- signature[signature$Signature == sig_name, "method"]

    geneSig <- switch(method,
        "GSVA" = geneSigGSVA(dat.icb = ICB_small_Mariathasan, sig = sig, sig.name = sig_name, missing.perc = 0.5, const.int = 0.001, n.cutoff = 15, sig.perc = 0.8, study = "ICB_small_Mariathasan"),
        "Weighted Mean" = geneSigMean(dat.icb = ICB_small_Mariathasan, sig = sig, sig.name = sig_name, missing.perc = 0.5, const.int = 0.001, n.cutoff = 15, sig.perc = 0.8, study = "ICB_small_Mariathasan"),
        "ssGSEA" = geneSigssGSEA(dat.icb = ICB_small_Mariathasan, sig = sig, sig.name = sig_name, missing.perc = 0.5, const.int = 0.001, n.cutoff = 15, sig.perc = 0.8, study = "ICB_small_Mariathasan"),
        "Specific Algorithm" = switch(sig_name,
            "COX-IS_Bonavita" = geneSigCOX_IS(dat.icb = ICB_small_Mariathasan, sig = sig, sig.name = sig_name, missing.perc = 0.5, const.int = 0.001, n.cutoff = 15, sig.perc = 0.8, study = "ICB_small_Mariathasan"),
            "IPS_Charoentong" = geneSigIPS(dat.icb = ICB_small_Mariathasan, sig = sig, sig.name = sig_name, missing.perc = 0.5, const.int = 0.001, n.cutoff = 15, sig.perc = 0.8, study = "ICB_small_Mariathasan"),
            "PredictIO_Bareche" = geneSigPredictIO(dat.icb = ICB_small_Mariathasan, sig = sig, sig.name = sig_name, missing.perc = 0.5, const.int = 0.001, n.cutoff = 15, sig.perc = 0.8, study = "ICB_small_Mariathasan"),
            "IPRES_Hugo" = geneSigIPRES(dat.icb = ICB_small_Mariathasan, sig = sig, sig.name = sig_name, missing.perc = 0.5, const.int = 0.001, n.cutoff = 15, sig.perc = 0.8, study = "ICB_small_Mariathasan"),
            "PassON_Du" = geneSigPassON(dat.icb = ICB_small_Mariathasan, sig = sig, sig.name = sig_name, missing.perc = 0.5, const.int = 0.001, n.cutoff = 15, sig.perc = 0.8, study = "ICB_small_Mariathasan"),
            "IPSOV_Shen" = geneSigIPSOV(dat.icb = ICB_small_Mariathasan, sig = sig, sig.name = sig_name, missing.perc = 0.5, const.int = 0.001, n.cutoff = 15, sig.perc = 0.8, study = "ICB_small_Mariathasan")
        )
    )

    if (sum(!is.na(geneSig)) == 0) geneSig <- rep(NA, ncol(ICB_small_Mariathasan))
    geneSig
})

geneSig.score <- do.call(rbind, geneSig.score)
rownames(geneSig.score) <- substr(basename(GeneSig_list), 1, nchar(basename(GeneSig_list)) - 4)
remove <- which(is.na(rowSums(geneSig.score)))
if (length(remove) > 0) geneSig.score <- geneSig.score[-remove, ]
write.csv(geneSig.score, file = "ICB_small_Mariathasan_GeneSigScore.csv", row.names = FALSE)
