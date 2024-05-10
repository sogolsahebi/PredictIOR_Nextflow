#!/bin/bash -ue
Rscript -e '
library(SummarizedExperiment);
library(dplyr);
load("ICB_Wolf__Breast__IO+chemo.rda");  
se <- get(ls()[1]);  
expr <- assay(se);  
clin <- colData(se) %>% as.data.frame();  
write.csv(expr, "expr_data.csv", row.names = FALSE);  
write.csv(clin, "clin_data.csv", row.names = FALSE);  
'
