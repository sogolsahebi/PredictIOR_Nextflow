#!/usr/bin/env Rscript
library(SummarizedExperiment)
library(PredictioR)
library(survival)
library(meta)
library(GSVA)

list_rda <- lapply(list.files(path = 'data', pattern = '*.rda', full.names = TRUE), function(file) {
    load(file)
})

expr <- list('ICB_Liu' = ICB_small_Liu, 'ICB_Padron' = ICB_small_Padron, 'ICB_Hugo' = ICB_small_Hugo, 
          'ICB_Mariathasan' = ICB_small_Mariathasan, 'ICB_Nathanson' = ICB_small_Nathanson, 
          'ICB_Riaz' = ICB_small_Riaz, 'ICB_Miao' = ICB_small_Miao, 'ICB_Van_Allen' = ICB_small_Van_Allen)

cancer_type <- c('Melanoma', 'Pancreas', 'Melanoma', 'Bladder', 'Melanoma', 'Melanoma', 'Kidney', 'Melanoma')
treatment_type <- c('PD-1/PD-L1', 'PD-1/PD-L1', 'PD-1/PD-L1', 'PD-1/PD-L1', 'CTLA4', 'IO+combo', 'PD-1/PD-L1', 'CTLA4')

geneSig_os <- lapply(1:length(expr), function(k){

   geneSig <- geneSigGSVA(dat.icb = expr[[k]], 
                          sig = CYT_Rooney,
                          sig.name = 'CYT_Rooney',
                          missing.perc = 0.5,
                          const.int = 0.001, 
                          n.cutoff = 15,
                          sig.perc = 0.8, 
                          study = names(expr)[k])

   if(sum(!is.na(geneSig)) > 0){

     res <- geneSigSurvDicho(dat.icb = expr[[k]],
                             geneSig = geneSig[1,],
                             time.censor = 24,
                             n.cutoff = 15,
                             study =  names(expr)[k],
                             surv.outcome = "OS",
                             n0.cutoff = 5,
                             n1.cutoff = 5,
                             sig.name = 'CYT_Rooney',
                             method = 'median',
                             var.type = FALSE,
                             cancer.type = cancer_type[k],
                             treatment = treatment_type[k])

   } else {

     res <- data.frame( Outcome = "OS",
                        Gene = NA, 
                        Study = names(expr)[k],
                        Coef = NA,
                        SE = NA,
                        N = NA,
                        Pval = NA,
                        Cancer_type= NA,
                        Treatment = NA) 
   }

   rownames(res) <- NULL

   res

 })

geneSig_os <- do.call(rbind, geneSig_os)
geneSig_os$FDR <- p.adjust(geneSig_os$Pval, method = "BH")
geneSig_os <- geneSig_os[order(geneSig_os$Pval, decreasing = FALSE), ]


# Now Aggregating Associations(OS) through Meta-analysis (Pan-cancer)

res <- metafun(coef = geneSig_os$Coef, 
           se = geneSig_os$SE, 
           study  = geneSig_os$Study, 
           pval = geneSig_os$Pval, 
           n = geneSig_os$N, 
           cancer.type = geneSig_os$Cancer_type,
           treatment = geneSig_os$Treatment,
           feature = unique(geneSig_os$Gene),
           cancer.spec = FALSE,
           treatment.spec = FALSE)

write.csv(geneSig_os, file = "ICB_small_Mariathasan_GeneSigAssociationOS.csv", row.names = FALSE)
write.csv(res, file = "ICB_small_Mariathasan_signature_Associations(OS)_pancancer.csv", row.names = FALSE)
