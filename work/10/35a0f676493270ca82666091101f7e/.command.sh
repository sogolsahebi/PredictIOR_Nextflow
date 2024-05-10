#!/bin/bash -ue
Rscript -e "library(SummarizedExperiment); se <- readRDS('\"ICB_Ravi__Lung__PD-\(L\)1.rda\"'); expr <- assay(se, 'expr'); write.table(head(expr), file='\"expr_head_ICB_Ravi__Lung__PD-(L)1.txt\"', quote=FALSE, sep='\t')"
