#!/usr/bin/env Rscript
source('/R/load_libraries.R')

signature <- read.csv("signature_information.csv")
signature$Signature <- as.character(signature$signature)
signature$method <- as.character(signature$method)

dir_GeneSig <- './SIG_data'
GeneSig_list <- list.files(path = dir_GeneSig, pattern = '*.rda', full.names = TRUE)

load("ICB_small_Mariathasan.rda")

compute_signature_score <- function(method, dat_icb, sig, sig_name) {
    switch(method,
        "GSVA" = geneSigGSVA(dat.icb = dat_icb, sig = sig, sig.name = sig_name, missing.perc = 0.5, const.int = 0.001, n.cutoff = 15, sig.perc = 0.8, study = "ICB_small_Mariathasan"),
        "Weighted Mean" = geneSigMean(dat.icb = dat_icb, sig = sig, sig.name = sig_name, missing.perc = 0.5, const.int = 0.001, n.cutoff = 15, sig.perc = 0.8, study = "ICB_small_Mariathasan"),
        "ssGSEA" = geneSigssGSEA(dat.icb = dat_icb, sig = sig, sig.name = sig_name, missing.perc = 0.5, const.int = 0.001, n.cutoff = 15, sig.perc = 0.8, study = "ICB_small_Mariathasan"),
        geneSigSpecific(dat.icb = dat_icb, sig = sig, sig.name = sig_name, method = method, study = "ICB_small_Mariathasan") 
    )
}

geneSig.score <- lapply(GeneSig_list, function(gene_sig_file) {
    load(gene_sig_file)
    sig_name <- substr(basename(gene_sig_file), 1, nchar(basename(gene_sig_file))-4)
    method <- signature[signature$Signature == sig_name, "method"]

    geneSig <- compute_signature_score(method, ICB_small_Mariathasan, sig, sig_name)

    if (sum(!is.na(geneSig)) == 0) {
        geneSig <- rep(NA, ncol(ICB_small_Mariathasan))
    }
    geneSig
})

geneSig.score <- do.call(rbind, geneSig.score)
rownames(geneSig.score) <- substr(basename(GeneSig_list), 1, nchar(basename(GeneSig_list)) - 4)
remove <- which(is.na(rowSums(geneSig.score)))
if (length(remove) > 0) geneSig.score <- geneSig.score[-remove, ]
write.csv(geneSig.score, file = "ICB_small_Mariathasan_GeneSigScore.csv", row.names = FALSE)
