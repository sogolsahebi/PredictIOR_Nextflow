#!/usr/bin/env Rscript
source('/R/load_libraries.R')

load("ICB_small_Padron.rda")
geneSig.score <- read.csv("ICB_small_Padron_GeneSigScore.csv")

print("Loaded geneSig.score:")
print(dim(geneSig.score))
print(head(geneSig.score))

res.logreg <- lapply(1:nrow(geneSig.score), function(k){
    sig_name <- rownames(geneSig.score)[k]
    geneSig_vector <- as.numeric(geneSig.score[k,])

    print(paste("Processing signature:", sig_name))
    print(paste("Length of geneSig_vector:", length(geneSig_vector)))
    print(paste("Number of NAs in geneSig_vector:", sum(is.na(geneSig_vector))))

    res <- geneSigLogReg(dat.icb = ICB_small_Padron,
                        geneSig = geneSig_vector,
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
