#!/usr/bin/env Rscript
source('/R/load_libraries.R')

# load signature_information.csv, extract the methods (like GSVA, WeightedMean or PredictIO)
signature <- read.csv("signature_information.csv")
signature$Signature <- as.character(signature$signature)
signature$method <- as.character(signature$method)

# extract all Signitures file

dir_GeneSig <- "./SIG_data"
GeneSig_files <- list.files(path = dir_GeneSig, pattern = '*.rda', full.names = TRUE)

GeneSig_list <- lapply(GeneSig_files, function(file_path) {
sig <- load(file_path)  # This assumes each .rda file contains one object, typically a data frame or list
get(sig)
})

# load the summarized experiments
load("ICB_small_Mariathasan.rda")

geneSig.score <- lapply(1:length(GeneSig_list), function(i) { 

    print(paste(i, GeneSig_list[i], sep="/"))
    load(file.path(dir_GeneSig, GeneSig_list[i]))
    sig_name <- substr(GeneSig_list[i], 1, nchar(GeneSig_list[i])-4)

    if (signature[signature$Signature == sig_name, "method"] == "GSVA") {
        geneSig <- geneSigGSVA(dat.icb = ICB_small_Mariathasan,
                               sig = sig,
                               sig.name = sig_name,
                               missing.perc = 0.5,
                               const.int = 0.001,
                               n.cutoff = 15,
                               sig.perc = 0.8,
                               study = "ICB_small_Mariathasan")

        if (sum(!is.na(geneSig)) > 0) {
            geneSig <- geneSig[1,]
        }     
    }

    if (signature[signature$Signature == sig_name, "method"] == "Weighted Mean") {
        geneSig <- geneSigMean(dat.icb = ICB_small_Mariathasan,
                               sig = sig,
                               sig.name = sig_name,
                               missing.perc = 0.5,
                               const.int = 0.001,
                               n.cutoff = 15,
                               sig.perc = 0.8,
                               study = "ICB_small_Mariathasan")
    }

    if (signature[signature$Signature == sig_name, "method"] == "ssGSEA") {
        geneSig <- geneSigssGSEA(dat.icb = ICB_small_Mariathasan,
                                 sig = sig,
                                 sig.name = sig_name,
                                 missing.perc = 0.5,
                                 const.int = 0.001,
                                 n.cutoff = 15,
                                 sig.perc = 0.8,
                                 study = "ICB_small_Mariathasan")

        if (sum(!is.na(geneSig)) > 0) {
            geneSig <- geneSig[1,]
        }   
    }

    if (signature[signature$Signature == sig_name, "method"] == "Specific Algorithm" & sig_name == "COX-IS_Bonavita") {
        geneSig <- geneSigCOX_IS(dat.icb = ICB_small_Mariathasan,
                                 sig = sig,
                                 sig.name = sig_name,
                                 missing.perc = 0.5,
                                 const.int = 0.001,
                                 n.cutoff = 15,
                                 sig.perc = 0.8,
                                 study = "ICB_small_Mariathasan")
    }

    if (signature[signature$Signature == sig_name, "method"] == "Specific Algorithm" & sig_name == "IPS_Charoentong") {
        geneSig <- geneSigIPS(dat.icb = ICB_small_Mariathasan,
                              sig = sig,
                              sig.name = sig_name,
                              missing.perc = 0.5,
                              const.int = 0.001,
                              n.cutoff = 15,
                              sig.perc = 0.8,
                              study = "ICB_small_Mariathasan")
    }

    if (signature[signature$Signature == sig_name, "method"] == "Specific Algorithm" & sig_name == "PredictIO_Bareche") {
        geneSig <- geneSigPredictIO(dat.icb = ICB_small_Mariathasan,
                                    sig = sig,
                                    sig.name = sig_name,
                                    missing.perc = 0.5,
                                    const.int = 0.001,
                                    n.cutoff = 15,
                                    sig.perc = 0.8,
                                    study = "ICB_small_Mariathasan")
    }

    if (signature[signature$Signature == sig_name, "method"] == "Specific Algorithm" & sig_name == "IPRES_Hugo") {
        geneSig <- geneSigIPRES(dat.icb = ICB_small_Mariathasan,
                                sig = sig,
                                sig.name = sig_name,
                                missing.perc = 0.5,
                                const.int = 0.001,
                                n.cutoff = 15,
                                sig.perc = 0.8,
                                study = "ICB_small_Mariathasan")
    }

    if (signature[signature$Signature == sig_name, "method"] == "Specific Algorithm" & sig_name == "PassON_Du") {
        geneSig <- geneSigPassON(dat.icb = ICB_small_Mariathasan,
                                 sig = sig,
                                 sig.name = sig_name,
                                 missing.perc = 0.5,
                                 const.int = 0.001,
                                 n.cutoff = 15,
                                 sig.perc = 0.8,
                                 study = "ICB_small_Mariathasan")
    }

    if (signature[signature$Signature == sig_name, "method"] == "Specific Algorithm" & sig_name == "IPSOV_Shen") {
        geneSig <- geneSigIPSOV(dat.icb = ICB_small_Mariathasan,
                                sig = sig,
                                sig.name = sig_name,
                                missing.perc = 0.5,
                                const.int = 0.001,
                                n.cutoff = 15,
                                sig.perc = 0.8,
                                study = "ICB_small_Mariathasan")
    }

    if (sum(!is.na(geneSig)) > 0) {
        return(geneSig)
    } else {
        return(rep(NA, ncol(expr)))
    }
})

geneSig.score <- do.call(rbind, geneSig.score)
rownames(geneSig.score) <- substr(GeneSig_list, 1, nchar(GeneSig_list) - 4)
remove <- which(is.na(rowSums(geneSig.score)))
if (length(remove) > 0) {
    geneSig.score <- geneSig.score[-remove, ]
}
write.csv(geneSig.score, file = "ICB_small_Mariathasan_GeneSigScore.csv", row.names = FALSE)
