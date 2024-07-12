#!/usr/bin/env Rscript
source('/R/load_libraries.R')

load("ICB_small_Van_Allen.rda")

geneSig.score <- read.csv("ICB_small_Van_Allen_GeneSigScore.csv")

res.all <- lapply(1:nrow(geneSig.score), function(k) {
    sig_name <- rownames(geneSig.score)[k]
    geneSig_vector <- as.numeric(geneSig.score[k, ])
    geneSig_vector <- geneSig_vector[!is.na(geneSig_vector)]

    res <- geneSigSurvCont(
        dat.icb = ICB_small_Van_Allen,
        geneSig = geneSig_vector,
        time.censor = 36,
        n.cutoff = 15,
        study = "ICB_small_Van_Allen",
        surv.outcome = "OS",
        sig.name = sig_name,
        cancer.type = "Pancreas",
        treatment = "PD-1/PD-L1"
    )

    res
})

res.all <- do.call(rbind, res.all)
res.all$FDR <- p.adjust(res.all$Pval, method="BH")
res.all <- res.all[order(res.all$FDR), ]

write.csv(res.all, file = "ICB_small_Van_Allen_os_GeneSig_association.csv", row.names = TRUE)
