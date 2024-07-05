#!/usr/bin/env Rscript
library(SummarizedExperiment)
library(PredictioR)
library(survival)
library(meta)
library(GSVA)
library(jsonlite) 

load("CYT_Rooney.rda") # this gives us 'sig'

# Load all .rda files 
lapply(list.files(path = 'ICB_data', pattern = '*.rda', full.names = TRUE), function(file) {
    load(file, envir = .GlobalEnv)
})

# Create a list of the loaded objects using their actual names
loaded_objects <- ls(pattern = "^ICB_small_")
expr <- mget(loaded_objects, envir = .GlobalEnv)

# Define the cancer types and treatment types vectors
cancer_type <- fromJSON('["Melanoma", "Pancreas", "Melanoma", "Bladder", "Melanoma", "Melanoma", "Kidney", "Melanoma"]')
treatment_type <- fromJSON('["PD-1/PD-L1", "PD-1/PD-L1", "PD-1/PD-L1", "PD-1/PD-L1", "CTLA4", "IO+combo", "PD-1/PD-L1", "CTLA4"]')

geneSig_os <- lapply(1:length(expr), function(k){

   geneSig <- geneSigGSVA(dat.icb = expr[[k]], 
                          sig = sig,
                          sig.name = "CYT_Rooney",
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
