#!/usr/bin/env Rscript
library(SummarizedExperiment)
library(PredictioR)
library(GSVA)

# Load ICB dataset
load("ICB_small_Mariathasan.rda")
expr <- assay(ICB_small_Mariathasan)
clin <- as.data.frame(colData(ICB_small_Mariathasan))

# List to store signature scores
signature_scores_list <- list()

# Process each signature method
for (sig_name in "[ADO_Sidders, APM_Thompson, APM_Wang, Bcell_Budczies, Bcell_Helmink, Blood_Friedlander, CCL5-CXCL9_Dangaj, CD39-CD8Tcell_Chow, CD8_Sade-Feldman, C-ECM_Chakravarthy, Chemokine_Messina, CIN25_Carter, CIN70_Carter, COX-IS_Bonavita, CRMA_Shukla, CYT_Davoli, CYT_Rooney, EMTstroma_Wang, EMT_Thompson, IFN_Ayers, ImmuneCells_Xiong, ImmuneScore_Roh, ImmunoScore_Hao, IMPRES_Auslander, IMS_Cui, Inflammatory_Thompson, IPRES_Hugo, IPS_Charoentong, IPSOV_Shen, IRG_Ayers, IRG_Yang, KDM5A_Wang, M1_Hwang, Macrophage_Danaher, MHC-II_Liu, MHC-I_Liu, MPS_Guijarro, Myeloid-DC_Helmink, NonResponse_Chen, PassON_Du, PDL1_Nishino, peri-Tcell_Hwang, PredictIO_Bareche, proliferation_Pabla, PTEN-MITF_Cabrita, Response_Chen, STAT1_Care, Tcell-exclusion_Arnon, Teff-IFNG_Fehrenbacher, Teff_McDermott, TGFB_Mariathasan, TIS_Damotte, TLS_Cabrita, TUMOR-VASCULATURE-UP_Lu, VIGex_Hernando-Calvo]") {
    # Load signature file
    sig_file <- file.path("SIG_data", paste0(sig_name, ".rda"))
    load(sig_file)
    sig <- get("sig", envir = .GlobalEnv)

    # Determine method for processing
    method <- switch(sig_name,
                     ADO_Sidders = 'GSVA',
                     APM_Thompson = 'Weighted Mean',
                     APM_Wang = 'GSVA',
                     default = 'GSVA'  # Default method if not specified
    )

    # Process signature based on method
    if (method == 'GSVA') {
        geneSig <- geneSigGSVA(dat.icb = expr,
                               sig = sig,
                               sig.name = sig_name,
                               missing.perc = 0.5,
                               const.int = 0.001,
                               n.cutoff = 15,
                               sig.perc = 0.8,
                               study = "ICB_small_Mariathasan")
    } else if (method == 'Weighted Mean') {
        geneSig <- geneSigMean(dat.icb = expr,
                               sig = sig,
                               sig.name = sig_name,
                               missing.perc = 0.5,
                               const.int = 0.001,
                               n.cutoff = 15,
                               sig.perc = 0.8,
                               study = "ICB_small_Mariathasan")
    }
    # Add more methods as needed

    # Store signature scores
    signature_scores_list[[sig_name]] <- geneSig
}

# Save all signature scores as CSV
write.csv(do.call(rbind, signature_scores_list), file = "./output/ICB_small_Mariathasan_SignatureScores.csv", row.names = TRUE)
