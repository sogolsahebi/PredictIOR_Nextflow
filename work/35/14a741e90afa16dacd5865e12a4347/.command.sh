#!/usr/bin/env Rscript
library(PredictIO)

expr <- read.csv("ICB_Ravi__Lung__PD-L1_expr.csv", row.names = 1)
clin <- read.csv("ICB_Ravi__Lung__PD-L1_clin.csv")

source('/R_scripts/getGeneAssociation.R')
source('/R_scripts/getHR.R')

cox_pfs <- lapply(1:100, function(i) {
  geneSurvCont(
    dat.icb = expr,
    clin = clin,
    time.censor = 36,
    missing.perc = 0.5,
    const.int = 0.001,
    n.cutoff = 15,
    feature = rownames(expr)[i],
    study = 'ICB_Ravi__Lung__PD-L1',
    surv.outcome = 'PFS',
    cancer.type = 'Lung',
    treatment = 'PD-L1'
  )
})

cox_pfs <- do.call(rbind, cox_pfs)

# Additional filtering and adjustments as needed
cox_pfs <- cox_pfs[!is.na(cox_pfs$Gene), ]
cox_pfs$FDR <- p.adjust(cox_pfs$Pval, method = "BH")

write.csv(cox_pfs, file = "ICB_Ravi__Lung__PD-L1_cox_pfs_results.csv", row.names = FALSE)
