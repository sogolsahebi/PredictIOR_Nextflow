#!/usr/bin/env Rscript
library(SummarizedExperiment)

expr <- read.csv("ICB_Ravi__Lung__PD-L1_expr.csv", row.names = 1)
clin <- read.csv("ICB_Ravi__Lung__PD-L1_clin.csv")

source('/R_scripts/getGeneAssociation.R')
source('/R_scripts/getHR.R')

results <- lapply(1:100, function(i) {
    geneSurvCont(
        dat.icb = expr,
        clin = clin,
        time.censor = 36,
        missing.perc = 0.5,
        const.int = 0.001,
        n.cutoff = 15,
        feature = rownames(expr)[i],
        study = 'ICB_Ravi__Lung__PD-L1',
        surv.outcome = 'OS',
        cancer.type = 'Lung',
        treatment = 'PD-L1'
    )
})

cox_os <- do.call(rbind, results)
# Additional filtering and adjustments as needed
cox_os <- cox_os[!is.na(cox_os$Gene), ]
cox_os$FDR <- p.adjust(cox_os$Pval, method = "BH")

write.csv(cox_os, file = "ICB_Ravi__Lung__PD-L1_cox_os_results.csv", row.names = FALSE)
