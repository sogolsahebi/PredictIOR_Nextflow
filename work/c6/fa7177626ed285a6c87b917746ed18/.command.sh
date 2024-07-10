#!/usr/bin/env Rscript
source('/R/load_libraries.R')

signature <- read.csv("signature_information.csv")
signature$Signature <- as.character(signature$signature)
signature$method <- as.character(signature$method)

print(head(signature[, c("Signature", "method")]))
print("methods")
print(signature$method)

print("Loaded signature information:")
print(head(signature))

dir_GeneSig <- './SIG_data'
GeneSig_list <- list.files(path = './SIG_data', pattern = '*.rda', full.names = TRUE)


print("Loaded dir_GeneSig")
print(length(GeneSig_list))


load("ICB_small_Mariathasan.rda")

geneSig.score <- lapply(1:length(GeneSig_list), function(i) {

    load(GeneSig_list[i])
    print("lest check sig type")
    print(class(sig))
    print(str(sig))

    sig_name <- substr(basename(GeneSig_list[i]), 1, nchar(basename(GeneSig_list[i]))-4)

    print(" see sig_name")
    print(sig_name)

    method <- signature[signature$Signature == sig_name, "method"]

    print(paste("Signature Name:", sig_name))
    print(paste("Method:", method))


    if (method == "GSVA") { geneSig <- geneSigGSVA(dat.icb = ICB_small_Mariathasan, 
                            sig = sig, 
                            sig.name = sig_name, 
                            missing.perc = 0.5, 
                            const.int = 0.001, 
                            n.cutoff = 15, 
                            sig.perc = 0.8)
    } else if (method == "Weighted Mean") {
        geneSig <- geneSigMean(dat.icb = ICB_small_Mariathasan, sig = sig, sig.name = sig_name, missing.perc = 0.5, const.int = 0.001, n.cutoff = 15, sig.perc = 0.8)
    } else if (method == "ssGSEA") {
        geneSig <- geneSigssGSEA(dat.icb = ICB_small_Mariathasan, sig = sig, sig.name = sig_name, missing.perc = 0.5, const.int = 0.001, n.cutoff = 15, sig.perc = 0.8)
    } else if (method == "Specific Algorithm") {
        geneSig <- switch(sig_name,
                          "COX-IS_Bonavita" = geneSigCOX_IS(dat.icb = ICB_small_Mariathasan, sig = sig, sig.name = sig_name, missing.perc = 0.5, const.int = 0.001, n.cutoff = 15, sig.perc = 0.8),
                          "IPS_Charoentong" = geneSigIPS(dat.icb = ICB_small_Mariathasan, sig = sig, sig.name = sig_name, missing.perc = 0.5, const.int = 0.001, n.cutoff = 15, sig.perc = 0.8),
                          geneSig)
    } else {
        print(paste("Unknown method for signature", sig_name, ":", method))
    }
    if (!is.null(geneSig) && sum(!is.na(geneSig)) > 0) {
        geneSig <- geneSig[1,]
    } else {
        geneSig <- rep(NA, ncol(ICB_small_Mariathasan))
    }
    return(geneSig)
})

geneSig.score <- do.call(rbind, geneSig.score)
rownames(geneSig.score) <- substr(GeneSig_files, 1, nchar(GeneSig_files) - 4)
remove <- which(is.na(rowSums(geneSig.score)))
if (length(remove) > 0) {
    geneSig.score <- geneSig.score[-remove, ]
}
write.csv(geneSig.score, file = "ICB_small_Mariathasan_GeneSigScore.csv", row.names = FALSE)
