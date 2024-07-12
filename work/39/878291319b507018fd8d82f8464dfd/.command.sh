#!/usr/bin/env Rscript
source('/R/load_libraries.R')

load("ICB_small_Liu.rda")

geneSig.score <- read.csv("ICB_small_Liu_GeneSigScore.csv", row.names = 1)

res.all <- lapply(1:nrow(geneSig.score), function(k) {
    sig_name <- rownames(geneSig.score)[k]
    geneSig_vector <- as.numeric(geneSig.score[k, ])
    geneSig_vector <- geneSig_vector[!is.na(geneSig_vector)]

    res <- geneSigSurvCont(
        dat.icb = ICB_small_Liu,
        geneSig = geneSig_vector,
        time.censor = 36,
        n.cutoff = 15,
        study = "ICB_small_Liu",
        surv.outcome = "OS",
        sig.name = sig_name,
        cancer.type = "Melanoma",
        treatment = "PD-1/PD-L1"
    )

    res
})

res.all <- do.call(rbind, res.all)
res.all$FDR <- p.adjust(res.all$Pval, method="BH")
res.all <- res.all[order(res.all$FDR), ]

write.csv(res.all, file = "ICB_small_Liu_os_GeneSig_association.csv", row.names = TRUE)
