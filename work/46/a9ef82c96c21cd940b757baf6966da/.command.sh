#!/usr/bin/env Rscript
source('/R/load_libraries.R')

geneSig.score <- read.csv("ICB_small_Mariathasan_GeneSigScore.csv")

# Ensure rownames are retained correctly after reading the CSV
rownames(geneSig.score) <- geneSig.score[,1]
geneSig.score <- geneSig.score[,-1]

res.all <- lapply(1:nrow(geneSig.score), function(k) {
    print(paste("Processing signature:", rownames(geneSig.score)[k]))

    geneSig <- geneSig.score[k, ]
    print(head(geneSig))

    # Check for NA values and handle them
    dat.icb <- get(dat.icb)  # Load your dataset here
    dat.icb <- dat.icb[!is.na(dat.icb$time) & !is.na(dat.icb$status), ]
    dat.icb$time <- as.numeric(as.character(dat.icb$time))
    dat.icb$status <- as.numeric(as.character(dat.icb$status))

    # Ensure no NA values are present after conversion
    if (any(is.na(dat.icb$time)) || any(is.na(dat.icb$status))) {
        stop("NA values found in time or status after conversion to numeric")
    }

    res <- geneSigSurvCont(
        dat.icb = dat.icb,
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
