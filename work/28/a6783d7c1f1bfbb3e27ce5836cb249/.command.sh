#!/usr/bin/env Rscript
library(SummarizedExperiment)
library(PredictioR)
library(GSVA)
library(jsonlite)

# Load icb rda file
load("ICB_small_Mariathasan.rda")
expr <- assay(ICB_small_Mariathasan)
clin <- as.data.frame(colData(ICB_small_Mariathasan))

# List all signature files
sig_files <- list.files(path = "SIG_data", pattern = "*.rda", full.names = TRUE)

# Define method mapping
method_mapping <- ${params.signature_methods as JSON}

signature_scores <- lapply(sig_files, function(sig_file) {
    load(sig_file) # Ensure this loads 'sig'
    sig_name <- sub(".rda", "", basename(sig_file))

    method <- method_mapping[[sig_name]]

    if (method == "GSVA") {
        geneSig <- geneSigGSVA(dat.icb = expr,
                               sig = sig,
                               sig.name = sig_name,
                               missing.perc = 0.5,
                               const.int = 0.001,
                               n.cutoff = 15,
                               sig.perc = 0.8,
                               study = "ICB_small_Mariathasan")
    } else if (method == "Weighted Mean") {
        geneSig <- geneSigMean(dat.icb = expr,
                               sig = sig,
                               sig.name = sig_name,
                               missing.perc = 0.5,
                               const.int = 0.001,
                               n.cutoff = 15,
                               sig.perc = 0.8,
                               study = "ICB_small_Mariathasan")
    } else if (method == "ssGSEA") {
        geneSig <- geneSigssGSEA(dat.icb = expr,
                                 sig = sig,
                                 sig.name = sig_name,
                                 missing.perc = 0.5,
                                 const.int = 0.001,
                                 n.cutoff = 15,
                                 sig.perc = 0.8,
                                 study = "ICB_small_Mariathasan")
    } else if (method == "Specific Algorithm") {
        if (sig_name == "COX-IS_Bonavita") {
            geneSig <- geneSigCOX_IS(dat.icb = expr,
                                     sig = sig,
                                     sig.name = sig_name,
                                     missing.perc = 0.5,
                                     const.int = 0.001,
                                     n.cutoff = 15,
                                     sig.perc = 0.8,
                                     study = "ICB_small_Mariathasan")
        } else if (sig_name == "IPS_Charoentong") {
            geneSig <- geneSigIPS(dat.icb = expr,
                                  sig = sig,
                                  sig.name = sig_name,
                                  missing.perc = 0.5,
                                  const.int = 0.001,
                                  n.cutoff = 15,
                                  sig.perc = 0.8,
                                  study = "ICB_small_Mariathasan")
        } else if (sig_name == "PredictIO_Bareche") {
            geneSig <- geneSigPredictIO(dat.icb = expr,
                                        sig = sig,
                                        sig.name = sig_name,
                                        missing.perc = 0.5,
                                        const.int = 0.001,
                                        n.cutoff = 15,
                                        sig.perc = 0.8,
                                        study = "ICB_small_Mariathasan")
        } else if (sig_name == "IPRES_Hugo") {
            geneSig <- geneSigIPRES(dat.icb = expr,
                                    sig = sig,
                                    sig.name = sig_name,
                                    missing.perc = 0.5,
                                    const.int = 0.001,
                                    n.cutoff = 15,
                                    sig.perc = 0.8,
                                    study = "ICB_small_Mariathasan")
        } else if (sig_name == "PassON_Du") {
            geneSig <- geneSigPassON(dat.icb = expr,
                                     sig = sig,
                                     sig.name = sig_name,
                                     missing.perc = 0.5,
                                     const.int = 0.001,
                                     n.cutoff = 15,
                                     sig.perc = 0.8,
                                     study = "ICB_small_Mariathasan")
        } else if (sig_name == "IPSOV_Shen") {
            geneSig <- geneSigIPSOV(dat.icb = expr,
                                    sig = sig,
                                    sig.name = sig_name,
                                    missing.perc = 0.5,
                                    const.int = 0.001,
                                    n.cutoff = 15,
                                    sig.perc = 0.8,
                                    study = "ICB_small_Mariathasan")
        }
    }

    if (sum(!is.na(geneSig)) > 0) {
        return(geneSig)
    } else {
        return(rep(NA, ncol(expr)))
    }
})

signature_scores <- do.call(rbind, signature_scores)
rownames(signature_scores) <- sub(".rda", "", basename(sig_files))
remove <- which(is.na(rowSums(signature_scores)))
if (length(remove) > 0) {
    signature_scores <- signature_scores[-remove, ]
}

save(signature_scores, file = "./output/ICB_small_Mariathasan_SignatureScores.rda")
write.csv(signature_scores, file = "./output/ICB_small_Mariathasan_SignatureScores.csv", row.names = TRUE)
