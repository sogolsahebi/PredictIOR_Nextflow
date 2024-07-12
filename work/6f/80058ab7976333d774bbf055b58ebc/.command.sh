#!/usr/bin/env Rscript
source('/R/load_libraries.R')

print("Reading signature information...")
signature <- read.csv("signature_information.csv")
print("Signature information loaded:")
print(head(signature))
print(paste("Number of rows in signature data frame:", nrow(signature)))

# Ensure the signature and method columns are properly read
if ("signature" %in% colnames(signature) && "method" %in% colnames(signature)) {
    signature$signature <- as.character(signature$signature)
    signature$method <- as.character(signature$method)
} else {
    stop("Required columns 'signature' and/or 'method' not found in the signature information.")
}

print("Structure of signature data frame:")
str(signature)

load("ICB_small_Padron.rda")

computeGeneSignature <- function(method, sig_name, sig) {
    if (!is.data.frame(sig)) {
        stop(paste("Error: sig is not a dataframe for signature", sig_name))
    }

    print(paste("Processing signature:", sig_name))
    print(paste("Using method:", method))

    if (method == "GSVA") {
        geneSig <- tryCatch({
            geneSig <- geneSigGSVA(dat.icb = ICB_small_Padron, sig = sig, sig.name = sig_name, missing.perc = 0.5, const.int = 0.001, n.cutoff = 15, sig.perc = 0.8, study = "ICB_small_Padron")
            if (sum(!is.na(geneSig)) > 0) {
                geneSig[1,]
            } else {
                rep(NA, ncol(ICB_small_Padron))
            }
        }, error = function(e) {
            print(paste("Error in GSVA for", sig_name, ":", e$message))
            return(rep(NA, ncol(ICB_small_Padron)))
        })
    } else if (method == "Weighted Mean") {
        geneSig <- tryCatch({
            geneSigMean(dat.icb = ICB_small_Padron, sig = sig, sig.name = sig_name, missing.perc = 0.5, const.int = 0.001, n.cutoff = 15, sig.perc = 0.8, study = "ICB_small_Padron")
        }, error = function(e) {
            print(paste("Error in Weighted Mean for", sig_name, ":", e$message))
            return(rep(NA, ncol(ICB_small_Padron)))
        })
    } else if (method == "ssGSEA") {
        geneSig <- tryCatch({
            geneSigssGSEA(dat.icb = ICB_small_Padron, sig = sig, sig.name = sig_name, missing.perc = 0.5, const.int = 0.001, n.cutoff = 15, sig.perc = 0.8, study = "ICB_small_Padron")
            if (sum(!is.na(geneSig)) > 0) {
                geneSig[1,]
            } else {
                rep(NA, ncol(ICB_small_Padron))
            }
        }, error = function(e) {
            print(paste("Error in ssGSEA for", sig_name, ":", e$message))
            return(rep(NA, ncol(ICB_small_Padron)))
        })
    } else if (method == "Specific Algorithm") {
        geneSig <- tryCatch({
            geneSigPassON(dat.icb = ICB_small_Padron, sig = sig, sig.name = sig_name, missing.perc = 0.5, const.int = 0.001, n.cutoff = 15, sig.perc = 0.8, study = "ICB_small_Padron")
        }, error = function(e) {
            print(paste("Error in Specific Algorithm for", sig_name, ":", e$message))
            return(rep(NA, ncol(ICB_small_Padron)))
        })
    } else {
        stop(paste("Unknown method:", method))
    }

    return(geneSig)
}

sig_files <- list.files(path = 'SIG_data', pattern = '*.rda', full.names = TRUE)

geneSig.score <- lapply(sig_files, function(file) {
    env <- new.env()
    load(file, envir = env)
    sig <- env$sig
    sig_name <- substr(basename(file), 1, nchar(basename(file)) - 4)
    method <- signature[signature$signature == sig_name, "method"]
    if (length(method) == 0) {
        print(paste("No method found for signature", sig_name))
        return(rep(NA, ncol(ICB_small_Padron)))
    }
    computeGeneSignature(method, sig_name, sig)
})

geneSig.score <- do.call(rbind, geneSig.score)
rownames(geneSig.score) <- substr(basename(sig_files), 1, nchar(basename(sig_files)) - 4)
remove <- which(is.na(rowSums(geneSig.score)))
if (length(remove) > 0) {
    geneSig.score <- geneSig.score[-remove, ]
}
write.csv(geneSig.score, file = "ICB_small_Padron_GeneSigScore.csv", row.names = FALSE)