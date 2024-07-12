#!/usr/bin/env Rscript
source('/R/load_libraries.R')

load("ICB_small_Padron.rda")

geneSig.score <- read.csv("ICB_small_Padron_GeneSigScore.csv")

print("Loaded geneSig.score:")
print(dim(geneSig.score))
print(head(geneSig.score))

res.all <- lapply(1:nrow(geneSig.score), function(k) {
    sig_name <- rownames(geneSig.score)[k]
    geneSig_vector <- as.numeric(geneSig.score[k, ])
    geneSig_vector <- as.numeric(geneSig.score[k, ])
    geneSig_vector <- geneSig_vector[!is.na(geneSig_vector)]

    print(paste("Processing signature:", sig_name))
    print(paste("Length of geneSig_vector:", length(geneSig_vector)))
    print(paste("Number of NAs in geneSig_vector:", sum(is.na(geneSig_vector))))

    res <- geneSigSurvCont(
        dat.icb = ICB_small_Padron,
        geneSig = geneSig_vector,
        time.censor = 36,
        n.cutoff = 15,
        study = "ICB_small_Padron",
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

write.csv(res.all, file = "ICB_small_Padron_os_GeneSig_association.csv", row.names = TRUE)
