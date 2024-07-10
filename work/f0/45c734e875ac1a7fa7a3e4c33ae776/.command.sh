#!/usr/bin/env Rscript
source('/R/load_libraries.R')

load("ICB_small_Mariathasan.rda")

geneSig.score <- read.csv("ICB_small_Mariathasan_GeneSigScore.csv")

# Ensure rownames are retained correctly after reading the CSV
rownames(geneSig.score) <- geneSig.score[,1]
geneSig.score <- geneSig.score[,-1]

res.all <- lapply(1:nrow(geneSig.score), function(k) {
    print(paste("Processing signature:", rownames(geneSig.score)[k]))

    geneSig <- geneSig.score[k, ]

    # Check for missing values in the data
    print(paste("Checking data for signature:", rownames(geneSig.score)[k]))
    missing_values <- sum(is.na(ICB_small_Mariathasan$time))
    print(paste("Number of missing values in time:", missing_values))
    missing_status <- sum(is.na(ICB_small_Mariathasan$status))
    print(paste("Number of missing values in time:", missing_status))

    # Handle missing values
    ICB_small_Mariathasan <- ICB_small_Mariathasan[!is.na(ICB_small_Mariathasan$time), ]

    res <- geneSigSurvCont(
        dat.icb = ICB_small_Mariathasan,
        geneSig = geneSig,
        time.censor = 24,
        n.cutoff = 15,
        study = "ICB_small_Mariathasan",
        surv.outcome = "OS",
        sig.name = rownames(geneSig.score)[k],
        cancer.type = "Bladder",
        treatment = "PD-1/PD-L1"
    )

    print(res)
    res
})

res.all <- do.call(rbind, res.all)
res.all$FDR <- p.adjust(res.all$Pval, method="BH")
res.all <- res.all[order(res.all$FDR), ]

write.csv(res.all, file = "ICB_small_Mariathasan_cox_OS_association.csv", row.names = FALSE)
