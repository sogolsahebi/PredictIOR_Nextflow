#!/usr/bin/env Rscript
source('/R/load_libraries.R')

geneSig.score <- read.csv("ICB_small_Mariathasan_GeneSigScore.csv")

res.all <- lapply(1:nrow(geneSig.score), function(k) {
    res <- geneSigSurvCont(
        dat.icb = ICB_small_Mariathasan,
        geneSig = geneSig.score[k, ],
        time.censor = 24,
        n.cutoff = 15,
        study = "ICB_small_Mariathasan",
        surv.outcome = "OS",
        sig.name = rownames(geneSig.score)[k],
        cancer.type = "Bladder",
        treatment = "PD-1/PD-L1"
    )
    res
})

res.all <- do.call(rbind, res.all)
res.all$FDR <- p.adjust(res.all$Pval, method="BH")
res.all <- res.all[order(res.all$FDR), ]

write.csv(res.all, file = "ICB_small_Mariathasan_cox_OS_association.csv", row.names = FALSE)
