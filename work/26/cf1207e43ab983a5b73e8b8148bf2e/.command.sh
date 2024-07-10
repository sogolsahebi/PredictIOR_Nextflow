#!/usr/bin/env Rscript
source('/R/load_libraries.R')

load("ICB_small_Mariathasan.rda")

geneSig.score <- read.csv("ICB_small_Mariathasan_GeneSigScore.csv", row.names = 1)

# Ensure rownames are retained correctly after reading the CSV
rownames(geneSig.score) <- geneSig.score[,1]
geneSig.score <- geneSig.score[,-1]

res.all <- lapply(1:nrow(geneSig.score), function(k) {
    print(paste("Processing signature:", rownames(geneSig.score)[k]))

    res <- geneSigSurvCont(
        dat.icb = ICB_small_Mariathasan,
        geneSig = as.numeric(geneSig.score[k, ]),
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
#res.all$FDR <- p.adjust(res.all$Pval, method="BH")
#res.all <- res.all[order(res.all$FDR), ]

write.csv(res.all, file = "ICB_small_Mariathasan_OS_association.csv", row.names = FALSE)