#!/usr/bin/env Rscript
source('/R/load_libraries.R')

# load signature_information.csv, extract the methods (like GSVA, WeightedMean or PredictIO)
signature <- read.csv("signature_information.csv")
signature$Signature <- as.character(signature$signature)
signature$method <- as.character(signature$method)

# extract all Signitures file

dir_GeneSig <- 'ADO_Sidders.rda APM_Thompson.rda APM_Wang.rda Bcell_Budczies.rda Bcell_Helmink.rda Blood_Friedlander.rda C-ECM_Chakravarthy.rda CCL5-CXCL9_Dangaj.rda CD39-CD8Tcell_Chow.rda CD8_Sade-Feldman.rda CIN25_Carter.rda CIN70_Carter.rda COX-IS_Bonavita.rda CRMA_Shukla.rda CYT_Davoli.rda CYT_Rooney.rda Chemokine_Messina.rda EMT_Thompson.rda EMTstroma_Wang.rda IFN_Ayers.rda IMPRES_Auslander.rda IMS_Cui.rda IPRES_Hugo.rda IPSOV_Shen.rda IPS_Charoentong.rda IRG_Ayers.rda IRG_Yang.rda ImmuneCells_Xiong.rda ImmuneScore_Roh.rda ImmunoScore_Hao.rda Inflammatory_Thompson.rda KDM5A_Wang.rda M1_Hwang.rda MHC-II_Liu.rda MHC-I_Liu.rda MPS_Guijarro.rda Macrophage_Danaher.rda Myeloid-DC_Helmink.rda NonResponse_Chen.rda PDL1_Nishino.rda PTEN-MITF_Cabrita.rda PassON_Du.rda PredictIO_Bareche.rda Response_Chen.rda STAT1_Care.rda TGFB_Mariathasan.rda TIS_Damotte.rda TLS_Cabrita.rda TUMOR-VASCULATURE-UP_Lu.rda Tcell-exclusion_Arnon.rda Teff-IFNG_Fehrenbacher.rda Teff_McDermott.rda VIGex_Hernando-Calvo.rda peri-Tcell_Hwang.rda proliferation_Pabla.rda'

GeneSig_files <- list.files(path = dir_GeneSig, pattern = '*.rda', full.names = TRUE)
sig_list <- lapply(GeneSig_files, function(file) {
    load(file)
    return(sig)
})

# load the summarized experiments
load("ICB_small_Mariathasan.rda")

geneSig.score <- lapply(1:length(GeneSig_list), function(i) { 

    print(paste(i, GeneSig_list[i], sep="/"))
    load(file.path(dir_GeneSig, GeneSig_list[i]))
    sig_name <- GeneSig_list[i]
    sig_name <- sig  

    if (signature[signature$Signature == sig_name, "method"] == "GSVA") {
        geneSig <- geneSigGSVA(dat.icb = ICB_small_Mariathasan,
                               sig = sig,
                               sig.name = sig_name,
                               missing.perc = 0.5,
                               const.int = 0.001,
                               n.cutoff = 15,
                               sig.perc = 0.8,
                               study = "ICB_small_Mariathasan")

        if (sum(!is.na(geneSig)) > 0) {
            geneSig <- geneSig[1,]
        }     
    }

    if (signature[signature$Signature == sig_name, "method"] == "Weighted Mean") {
        geneSig <- geneSigMean(dat.icb = ICB_small_Mariathasan,
                               sig = sig,
                               sig.name = sig_name,
                               missing.perc = 0.5,
                               const.int = 0.001,
                               n.cutoff = 15,
                               sig.perc = 0.8,
                               study = "ICB_small_Mariathasan")
    }

    if (signature[signature$Signature == sig_name, "method"] == "ssGSEA") {
        geneSig <- geneSigssGSEA(dat.icb = ICB_small_Mariathasan,
                                 sig = sig,
                                 sig.name = sig_name,
                                 missing.perc = 0.5,
                                 const.int = 0.001,
                                 n.cutoff = 15,
                                 sig.perc = 0.8,
                                 study = "ICB_small_Mariathasan")

        if (sum(!is.na(geneSig)) > 0) {
            geneSig <- geneSig[1,]
        }   
    }

    if (signature[signature$Signature == sig_name, "method"] == "Specific Algorithm" & sig_name == "COX-IS_Bonavita") {
        geneSig <- geneSigCOX_IS(dat.icb = ICB_small_Mariathasan,
                                 sig = sig,
                                 sig.name = sig_name,
                                 missing.perc = 0.5,
                                 const.int = 0.001,
                                 n.cutoff = 15,
                                 sig.perc = 0.8,
                                 study = "ICB_small_Mariathasan")
    }

    if (signature[signature$Signature == sig_name, "method"] == "Specific Algorithm" & sig_name == "IPS_Charoentong") {
        geneSig <- geneSigIPS(dat.icb = ICB_small_Mariathasan,
                              sig = sig,
                              sig.name = sig_name,
                              missing.perc = 0.5,
                              const.int = 0.001,
                              n.cutoff = 15,
                              sig.perc = 0.8,
                              study = "ICB_small_Mariathasan")
    }

    if (signature[signature$Signature == sig_name, "method"] == "Specific Algorithm" & sig_name == "PredictIO_Bareche") {
        geneSig <- geneSigPredictIO(dat.icb = ICB_small_Mariathasan,
                                    sig = sig,
                                    sig.name = sig_name,
                                    missing.perc = 0.5,
                                    const.int = 0.001,
                                    n.cutoff = 15,
                                    sig.perc = 0.8,
                                    study = "ICB_small_Mariathasan")
    }

    if (signature[signature$Signature == sig_name, "method"] == "Specific Algorithm" & sig_name == "IPRES_Hugo") {
        geneSig <- geneSigIPRES(dat.icb = ICB_small_Mariathasan,
                                sig = sig,
                                sig.name = sig_name,
                                missing.perc = 0.5,
                                const.int = 0.001,
                                n.cutoff = 15,
                                sig.perc = 0.8,
                                study = "ICB_small_Mariathasan")
    }

    if (signature[signature$Signature == sig_name, "method"] == "Specific Algorithm" & sig_name == "PassON_Du") {
        geneSig <- geneSigPassON(dat.icb = ICB_small_Mariathasan,
                                 sig = sig,
                                 sig.name = sig_name,
                                 missing.perc = 0.5,
                                 const.int = 0.001,
                                 n.cutoff = 15,
                                 sig.perc = 0.8,
                                 study = "ICB_small_Mariathasan")
    }

    if (signature[signature$Signature == sig_name, "method"] == "Specific Algorithm" & sig_name == "IPSOV_Shen") {
        geneSig <- geneSigIPSOV(dat.icb = ICB_small_Mariathasan,
                                sig = sig,
                                sig.name = sig_name,
                                missing.perc = 0.5,
                                const.int = 0.001,
                                n.cutoff = 15,
                                sig.perc = 0.8,
                                study = "ICB_small_Mariathasan")
    }

    if (sum(!is.na(geneSig)) > 0) {
        return(geneSig)
    } else {
        return(rep(NA, ncol(expr)))
    }
})

geneSig.score <- do.call(rbind, geneSig.score)
rownames(geneSig.score) <- substr(GeneSig_list, 1, nchar(GeneSig_list) - 4)
remove <- which(is.na(rowSums(geneSig.score)))
if (length(remove) > 0) {
    geneSig.score <- geneSig.score[-remove, ]
}
write.csv(geneSig.score, file = "ICB_small_Mariathasan_GeneSigScore.csv", row.names = FALSE)
