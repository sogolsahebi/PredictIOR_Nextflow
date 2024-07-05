#!/usr/bin/env Rscript
library(SummarizedExperiment)
library(PredictioR)
library(GSVA)
library(jsonlite)

# Load icb rda file
load("ICB_small_Mariathasan.rda")
expr <- assay(ICB_small_Mariathasan)
clin <- as.data.frame(colData(ICB_small_Mariathasan))

# List all signature files
sig_files <- list.files(path = "SIG_data", pattern = "*.rda", full.names = TRUE)

# Define method mapping
method_mapping <- fromJSON('[ADO_Sidders:GSVA, APM_Thompson:Weighted Mean, APM_Wang:GSVA, Bcell_Budczies:GSVA, Bcell_Helmink:GSVA, Blood_Friedlander:GSVA, CCL5-CXCL9_Dangaj:GSVA, CD39-CD8Tcell_Chow:GSVA, CD8_Sade-Feldman:GSVA, C-ECM_Chakravarthy:GSVA, Chemokine_Messina:GSVA, CIN25_Carter:GSVA, CIN70_Carter:GSVA, COX-IS_Bonavita:Specific Algorithm, CRMA_Shukla:GSVA, CYT_Davoli:GSVA, CYT_Rooney:GSVA, EMTstroma_Wang:GSVA, EMT_Thompson:Weighted Mean, IFN_Ayers:GSVA, ImmuneCells_Xiong:GSVA, ImmuneScore_Roh:GSVA, ImmunoScore_Hao:GSVA, IMPRES_Auslander:GSVA, IMS_Cui:GSVA, Inflammatory_Thompson:GSVA, IPRES_Hugo:Specific Algorithm, IPS_Charoentong:Specific Algorithm, IPSOV_Shen:Specific Algorithm, IRG_Ayers:GSVA, IRG_Yang:GSVA, KDM5A_Wang:GSVA, M1_Hwang:GSVA, Macrophage_Danaher:GSVA, MHC-II_Liu:GSVA, MHC-I_Liu:GSVA, MPS_Guijarro:GSVA, Myeloid-DC_Helmink:GSVA, NonResponse_Chen:GSVA, PassON_Du:Specific Algorithm, PDL1_Nishino:GSVA, peri-Tcell_Hwang:GSVA, PredictIO_Bareche:Specific Algorithm, proliferation_Pabla:GSVA, PTEN-MITF_Cabrita:GSVA, Response_Chen:GSVA, STAT1_Care:GSVA, Tcell-exclusion_Arnon:GSVA, Teff-IFNG_Fehrenbacher:GSVA, Teff_McDermott:GSVA, TGFB_Mariathasan:GSVA, TIS_Damotte:GSVA, TLS_Cabrita:GSVA, TUMOR-VASCULATURE-UP_Lu:GSVA, VIGex_Hernando-Calvo:GSVA]')

signature_scores <- lapply(sig_files, function(sig_file) {
    load(sig_file) # Ensure this loads 'sig'
    sig_name <- sub(".rda", "", basename(sig_file))

    method <- method_mapping[[sig_name]]

    if (method == "GSVA") {
        geneSig <- geneSigGSVA(dat.icb = expr,
                               sig = sig,
                               sig.name = sig_name,
                               missing.perc = 0.5,
                               const.int = 0.001,
                               n.cutoff = 15,
                               sig.perc = 0.8,
                               study = "ICB_small_Mariathasan")
    } else if (method == "Weighted Mean") {
        geneSig <- geneSigMean(dat.icb = expr,
                               sig = sig,
                               sig.name = sig_name,
                               missing.perc = 0.5,
                               const.int = 0.001,
                               n.cutoff = 15,
                               sig.perc = 0.8,
                               study = "ICB_small_Mariathasan")
    } else if (method == "ssGSEA") {
        geneSig <- geneSigssGSEA(dat.icb = expr,
                                 sig = sig,
                                 sig.name = sig_name,
                                 missing.perc = 0.5,
                                 const.int = 0.001,
                                 n.cutoff = 15,
                                 sig.perc = 0.8,
                                 study = "ICB_small_Mariathasan")
    } else if (method == "Specific Algorithm") {
        if (sig_name == "COX-IS_Bonavita") {
            geneSig <- geneSigCOX_IS(dat.icb = expr,
                                     sig = sig,
                                     sig.name = sig_name,
                                     missing.perc = 0.5,
                                     const.int = 0.001,
                                     n.cutoff = 15,
                                     sig.perc = 0.8,
                                     study = "ICB_small_Mariathasan")
        } else if (sig_name == "IPS_Charoentong") {
            geneSig <- geneSigIPS(dat.icb = expr,
                                  sig = sig,
                                  sig.name = sig_name,
                                  missing.perc = 0.5,
                                  const.int = 0.001,
                                  n.cutoff = 15,
                                  sig.perc = 0.8,
                                  study = "ICB_small_Mariathasan")
        } else if (sig_name == "PredictIO_Bareche") {
            geneSig <- geneSigPredictIO(dat.icb = expr,
                                        sig = sig,
                                        sig.name = sig_name,
                                        missing.perc = 0.5,
                                        const.int = 0.001,
                                        n.cutoff = 15,
                                        sig.perc = 0.8,
                                        study = "ICB_small_Mariathasan")
        } else if (sig_name == "IPRES_Hugo") {
            geneSig <- geneSigIPRES(dat.icb = expr,
                                    sig = sig,
                                    sig.name = sig_name,
                                    missing.perc = 0.5,
                                    const.int = 0.001,
                                    n.cutoff = 15,
                                    sig.perc = 0.8,
                                    study = "ICB_small_Mariathasan")
        } else if (sig_name == "PassON_Du") {
            geneSig <- geneSigPassON(dat.icb = expr,
                                     sig = sig,
                                     sig.name = sig_name,
                                     missing.perc = 0.5,
                                     const.int = 0.001,
                                     n.cutoff = 15,
                                     sig.perc = 0.8,
                                     study = "ICB_small_Mariathasan")
        } else if (sig_name == "IPSOV_Shen") {
            geneSig <- geneSigIPSOV(dat.icb = expr,
                                    sig = sig,
                                    sig.name = sig_name,
                                    missing.perc = 0.5,
                                    const.int = 0.001,
                                    n.cutoff = 15,
                                    sig.perc = 0.8,
                                    study = "ICB_small_Mariathasan")
        }
    }

    if (sum(!is.na(geneSig)) > 0) {
        return(geneSig)
    } else {
        return(rep(NA, ncol(expr)))
    }
})

signature_scores <- do.call(rbind, signature_scores)
rownames(signature_scores) <- sub(".rda", "", basename(sig_files))
remove <- which(is.na(rowSums(signature_scores)))
if (length(remove) > 0) {
    signature_scores <- signature_scores[-remove, ]
}

save(signature_scores, file = "./output/ICB_small_Mariathasan_SignatureScores.rda")
write.csv(signature_scores, file = "./output/ICB_small_Mariathasan_SignatureScores.csv", row.names = TRUE)
