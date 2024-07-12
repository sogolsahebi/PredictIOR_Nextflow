#!/usr/bin/env Rscript
source('/R/load_libraries.R')

load("ICB_small_Padron.rda")
geneSig.score <- read.csv("ICB_small_Padron_GeneSigScore.csv")
geneSig.score <- geneSig.score[,-1]


res.logreg <- lapply(1:nrow(geneSig.score), function(k){
sig_name <- rownames(geneSig.score)[k]
res <- geneSigLogReg(dat.icb = ICB_small_Padron,
                    geneSig = as.numeric(geneSig.score[k,]),
                    n.cutoff = 10,
                    study =  "ICB_small_Padron",
                    sig.name = sig_name,
                    n0.cutoff = 3, 
                    n1.cutoff = 3,
                    cancer.type = "Pancreas",
                    treatment = "PD-1/PD-L1")

res
})

res.logreg <- do.call(rbind, res.logreg)
res.logreg$FDR <- p.adjust(res.logreg$Pval, method="BH")
# Save as CSV file
write.csv(res.logreg, file = "ICB_small_Padron_GeneSig_Response.csv", row.names = TRUE)
