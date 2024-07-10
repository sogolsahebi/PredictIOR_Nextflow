#!/usr/bin/env Rscript
source('/R/load_libraries.R')

# load signature_information.csv, extract the methods (like GSVA, WeightedMean or PredictIO)
signature <- read.csv("signature_information.csv")
signature$Signature <- as.character(signature$signature)
signature$method <- as.character(signature$method)
print("Loaded signature information:")
print(head(signature))
print(signature$method)

# extract all Signatures file".rda"
dir_GeneSig <- "./SIG_data"
GeneSig_files <- list.files(dir_GeneSig)[order(list.files(dir_GeneSig))]

print("Loaded Genesig:")
print(length(GeneSig_files ))

# Assuming each .rda file contains one object, load and return the object by name.
GeneSig_list <- lapply(GeneSig_files, function(file_path) {
    load(file_path)
    data_name <- gsub(".rda", "", basename(file_path))
    get(data_name)
})

# load the summarized experiments
load("ICB_small_Mariathasan.rda")

# compute signatures scores based on different methods specified in signature_information
geneSig.score <- lapply(1:length(GeneSig_list), function(i) { 
    sig_name <- substr(GeneSig_list[i], 1, nchar(GeneSig_list[i])-4)
    method <- signature[signature$Signature == sig_name, "method"]
    geneSig <- NULL

    print(paste("Signature Name:", sig_name))
    print(paste("Method:", method))

    if (method == "GSVA") {
        geneSig <- geneSigGSVA(dat.icb = ICB_small_Mariathasan, sig = GeneSig_list[[i]], sig.name = sig_name, missing.perc = 0.5, const.int = 0.001, n.cutoff = 15, sig.perc = 0.8)
    } else if (method == "Weighted Mean") {
        geneSig <- geneSigMean(dat.icb = ICB_small_Mariathasan, sig = GeneSig_list[[i]], sig.name = sig_name, missing.perc = 0.5, const.int = 0.001, n.cutoff = 15, sig.perc = 0.8)
    } else if (method == "ssGSEA") {
        geneSig <- geneSigssGSEA(dat.icb = ICB_small_Mariathasan, sig = GeneSig_list[[i]], sig.name = sig_name, missing.perc = 0.5, const.int = 0.001, n.cutoff = 15, sig.perc = 0.8)
    } else if (method == "Specific Algorithm") {
        # Implement specific algorithm based on sig_name if applicable
        geneSig <- switch(sig_name,
                          "COX-IS_Bonavita" = geneSigCOX_IS(dat.icb = ICB_small_Mariathasan, sig = GeneSig_list[[i]], sig.name = sig_name, missing.perc = 0.5, const.int = 0.001, n.cutoff = 15, sig.perc = 0.8),
                          "IPS_Charoentong" = geneSigIPS(dat.icb = ICB_small_Mariathasan, sig = GeneSig_list[[i]], sig.name = sig_name, missing.perc = 0.5, const.int = 0.001, n.cutoff = 15, sig.perc = 0.8),
                          # add other specific algorithms here
                          geneSig)
    }
    if (!is.null(geneSig) && sum(!is.na(geneSig)) > 0) {
        geneSig <- geneSig[1,]
    } else {
        geneSig <- rep(NA, ncol(expr))
    }
    return(geneSig)
})

geneSig.score <- do.call(rbind, geneSig.score)
rownames(geneSig.score) <- substr(GeneSig_list, 1, nchar(GeneSig_list) - 4)
remove <- which(is.na(rowSums(geneSig.score)))
if (length(remove) > 0) {
    geneSig.score <- geneSig.score[-remove, ]
}
write.csv(geneSig.score, file = "ICB_small_Mariathasan_GeneSigScore.csv", row.names = FALSE)
