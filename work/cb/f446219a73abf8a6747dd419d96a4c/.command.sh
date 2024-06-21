#!/usr/bin/env Rscript
library(SummarizedExperiment)
library(PredictioR)

expr <- read.csv("ICB_small_Mariathasan_expr.csv", row.names = 1)
clin <- read.csv("ICB_small_Mariathasan_clin.csv")
cancer_type <- names( table(clin$cancer_type )[ table(clin$cancer_type ) >= 15 ] )


survCont( status = clin$event_occurred_os[ clin$cancer_type %in% cancer_type ] ,
      time = clin$survival_time_os[ clin$cancer_type %in% cancer_type ] ,
      time.censor = 36 , 
      var = as.numeric( scale( expr["CXCL9" , ] )))

#cox_os <- do.call(rbind, cox_os)

# Additional filtering 
# cox_os <- cox_os[!is.na(cox_os$Gene), ]
# cox_os$FDR <- p.adjust(cox_os$Pval, method = "BH")

#write.csv(cox_os, file = "ICB_small_Mariathasan_cox_os_genes.csv", row.names = FALSE)
