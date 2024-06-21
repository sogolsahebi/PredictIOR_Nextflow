#!/usr/bin/env Rscript
library(SummarizedExperiment)
library(PredictioR)
library(survival)

expr <- read.csv("ICB_small_Mariathasan_expr.csv", row.names = 1)
clin <- read.csv("ICB_small_Mariathasan_clin.csv")
cancer_type <- names(table(clin$cancer_type)[table(clin$cancer_type) >= 15])

# load "PredictIO_Bareche.rda"
load("PredictIO_Bareche.rda")

genes <- PredictIO_Bareche$gene_name
cox_os <- geneSurvCont(dat.icb = ICB_small_Mariathasan,
                    time.censor = 36,
                    missing.perc = 0.5,
                    const.int = 0.001,
                    n.cutoff = 15,
                    feature = genes,
                    study = ICB_small_Mariathasan,
                    surv.outcome = 'OS',
                    cancer.type = cancer_type,
                    treatment = 'PD-1/PD-L1')

# Additional filtering 
cox_os <- cox_os[order(cox_os$FDR <- p.adjust(cox_os$Pval, method = "BH")), ]
write.csv(cox_os, file = "ICB_small_Mariathasan_cox_os_barechegenes.csv", row.names = FALSE)
