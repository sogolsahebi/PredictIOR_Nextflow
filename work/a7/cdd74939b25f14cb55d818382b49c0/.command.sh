#!/usr/bin/env Rscript
source('/R/load_libraries.R')

load("ICB_small_Van_Allen.rda")
geneSig.score <- read.csv("ICB_small_Van_Allen_GeneSigScore.csv")


res.logreg <- lapply(1:nrow(geneSig.score), function(k){
sig_name <- rownames(geneSig.score)[k]
geneSig_vector <- as.numeric(geneSig.score[k, ])
geneSig_vector <- geneSig_vector[!is.na(geneSig_vector)]

res <- geneSigLogReg(dat.icb = ICB_small_Van_Allen,
                    geneSig = geneSig_vector,
                    n.cutoff = 10,
                    study =  "ICB_small_Van_Allen",
                    sig.name = sig_name,
                    n0.cutoff = 3, 
                    n1.cutoff = 3,
                    cancer.type = "Melanoma",
                    treatment = "PD-1/PD-L1")

res
})

res.logreg <- do.call(rbind, res.logreg)
res.logreg$FDR <- p.adjust(res.logreg$Pval, method="BH")
# Save as CSV file
write.csv(res.logreg, file = "ICB_small_Van_Allen_GeneSig_Response.csv", row.names = TRUE)