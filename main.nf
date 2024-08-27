#!/usr/bin/env nextflow
nextflow.enable.dsl=2

params.icb_data_dir = './ICB_data'
params.sig_data_dir = './SIG_data'
params.sig_summery_dir = './sig_summery_info' 
params.out_dir = './output/main_output'
params.gene_name = "CXCL9"

// Define cancer type and treatment for each study
cancer_type_map = [
    'ICB_small_Hugo'        : 'Melanoma',
    'ICB_small_Liu'         : 'Melanoma',
    'ICB_small_Miao'        : 'Kidney',
    'ICB_small_Nathanson'   : 'Melanoma',
    'ICB_small_Padron'      : 'Pancreas',
    'ICB_small_Riaz'        : 'Melanoma',
    'ICB_small_Van_Allen'   : 'Melanoma',
    'ICB_small_Mariathasan' : 'Bladder'
]

treatment_map = [
    'ICB_small_Hugo'        : 'PD-1/PD-L1',
    'ICB_small_Liu'         : 'PD-1/PD-L1',
    'ICB_small_Miao'        : 'PD-1/PD-L1',
    'ICB_small_Nathanson'   : 'CTLA4',
    'ICB_small_Padron'      : 'PD-1/PD-L1',
    'ICB_small_Riaz'        : 'PD-1/PD-L1',
    'ICB_small_Van_Allen'   : 'CTLA4',
    'ICB_small_Mariathasan' : 'PD-1/PD-L1'
]

/*
Note:
Define the cancer type and treatment for each study.
Another input example would be:

cancer_type_map = [
    'ICB_small_Liu' : 'Melanoma'
]

treatment_map = [
    'ICB_small_Liu' : 'PD-1/PD-L1'
]
*/

log.info """
P R E D I C T I O - N F   P I P E L I N E (Gene Level and Signature Level Analysis)
================================================================================
ICB Data Directory        : ${params.icb_data_dir}
Signature Data Directory  : ${params.sig_data_dir}
Output Directory          : ${params.out_dir}
""".stripIndent()

/*
========================================================
SECTION: Load Immunotherapy Datasets
========================================================
*/

/*
Public clinical datasets from GitHub or ORCESTRA. RNA profiles are log2-transformed TPM data from protein-coding genes, filtering out genes with zero expression in at least 50% of samples. Only studies with at least 20 patients are included.

For example, the Padron dataset includes RNA expression, clinical data, and gene metadata for 45 patients with 18,459 protein-coding genes, focused on Pancreas cancer and PD-1/PD-L1 treatment.

Links:
- GitHub: https://github.com/bhklab/PredictioR/tree/main/data
- ORCESTRA: https://www.orcestra.ca/clinical_icb
*/

// Load RDA data and extract expression, clinical data, and annotation data
process LoadAndExtractData {
    tag "${study_id}"
    container 'bhklab/nextflow-env:latest'
    publishDir "${params.out_dir}/${study_id}", mode: 'copy'

    input:
    tuple val(study_id), path(rda_file)

    output:
    tuple val(study_id), path("${study_id}_expr.csv"), path("${study_id}_clin.csv"), path("${study_id}_annot.csv")

    script:
    """
    #!/usr/bin/env Rscript
    source('/R/load_libraries.R')

    load("${rda_file}")

    # Extract expression data
    expr <- assay(${study_id})
    # Extract clinical data
    clin <- as.data.frame(colData(${study_id}))
    # Extract annotation data
    annot <- as.data.frame(rowData(${study_id}))
    
    # Write data to CSV files
    write.csv(expr, "${study_id}_expr.csv", row.names = TRUE)
    write.csv(clin, "${study_id}_clin.csv", row.names = FALSE)
    write.csv(annot, "${study_id}_annot.csv", row.names = TRUE)
    """
}

// Notes:
// For gene-level analysis in the R script, use dat.icb in two ways:
// 1. dat.icb = expr (data frame) with clin = clin (data frame) for clinical data.
// 2. Load the RDA file with load("${rda_file}"), then use dat.icb = '${study_id}' (SummarizedExperiment object) with clin = NULL.

/*
========================================================
SECTION: Biomarkers and Immunotherapy Response Association 
========================================================
*/

// Assessing the association of specific biomarkers with immunotherapy response (R vs NR) and survival (OS and PFS). P-values are corrected using the Benjamini-Hochberg (FDR) method, with significance set at p-values or FDR ≤ 5%.

// Gene association analysis for OS
process GeneAssociationOS {
    tag "${study_id}"
    container 'bhklab/nextflow-env:latest'
    publishDir "${params.out_dir}/${study_id}", mode: 'copy'

    input:
    tuple val(study_id), path(expr_file), path(clin_file), val(cancer_type), val(treatment), val(genes)

    output:
    path("${study_id}_cox_os.csv")

    script:
    """
    #!/usr/bin/env Rscript
    source('/R/load_libraries.R')
    
    # Read expression and clinical data from CSV files
    expr <- read.csv("${expr_file}", row.names = 1)
    clin <- read.csv("${clin_file}")

    # Perform gene association analysis for OS
    cox_result <- geneSurvCont(
        dat.icb = expr,
        clin = clin,
        time.censor = 36,
        missing.perc = 0.5,
        const.int = 0.001,
        n.cutoff = 15,
        feature = ${genes},
        study = "${study_id}",
        surv.outcome = "OS",
        cancer.type = "${cancer_type}",
        treatment = "${treatment}"
    )

    # Adjust p-values for multiple testing
    cox_result\$FDR <- p.adjust(cox_result\$Pval, method = "BH")
    cox_result <- cox_result[order(cox_result\$FDR), ]
    # Write results to CSV file
    write.csv(cox_result, file = "${study_id}_cox_os.csv", row.names = FALSE)
    """
}

// Gene association analysis for PFS
process GeneAssociationPFS {
    tag "${study_id}"
    container 'bhklab/nextflow-env:latest'
    publishDir "${params.out_dir}/${study_id}", mode: 'copy'

    input:
    tuple val(study_id), path(expr_file), path(clin_file), val(cancer_type), val(treatment), val(genes)

    output:
    path("${study_id}_cox_pfs.csv")

    script:
    """
    #!/usr/bin/env Rscript
    source('/R/load_libraries.R')
    
    # Read expression and clinical data from CSV files
    expr <- read.csv("${expr_file}", row.names = 1)
    clin <- read.csv("${clin_file}")

    # Perform gene association analysis for PFS
    cox_result <- geneSurvCont(
        dat.icb = expr,
        clin = clin,
        time.censor = 24,
        missing.perc = 0.5,
        const.int = 0.001,
        n.cutoff = 15,
        feature = ${genes},
        study = "${study_id}",
        surv.outcome = "PFS",
        cancer.type = "${cancer_type}",
        treatment = "${treatment}"
    )

    # Adjust p-values for multiple testing
    cox_result\$FDR <- p.adjust(cox_result\$Pval, method = "BH")
    cox_result <- cox_result[order(cox_result\$FDR), ]
    # Write results to CSV file
    write.csv(cox_result, file = "${study_id}_cox_pfs.csv", row.names = FALSE)
    """
}

// Gene association analysis for response
process GeneAssociationResponse {
    tag "${study_id}"
    container 'bhklab/nextflow-env:latest'
    publishDir "${params.out_dir}/${study_id}", mode: 'copy'

    input:
    tuple val(study_id), path(expr_file), path(clin_file), val(cancer_type), val(treatment), val(genes)

    output:
    path("${study_id}_logregResponse.csv")

    script:
    """
    #!/usr/bin/env Rscript
    source('/R/load_libraries.R')

    # Read expression and clinical data from CSV files
    expr <- read.csv("${expr_file}", row.names = 1)
    clin <- read.csv("${clin_file}")

    # Perform gene association analysis for response
    logreg <- geneLogReg(
        dat.icb = expr,
        clin = clin,
        missing.perc = 0.5,
        const.int = 0.001,
        n.cutoff = 15,
        feature = ${genes},
        study = "${study_id}", 
        n0.cutoff = 3, 
        n1.cutoff = 3,
        cancer.type = "${cancer_type}",
        treatment = "${treatment}"
    )

    # Adjust P-values and sort by FDR
    logreg\$FDR <- p.adjust(logreg\$Pval, method = "BH")
    logreg <- logreg[order(logreg\$FDR), ]
    
    # Save as CSV file
    write.csv(logreg, file = "${study_id}_logregResponse.csv", row.names = FALSE)
    """
}

/*
========================================================
SECTION: Signature Level Analysis
========================================================
*/

// Evaluation of over 50 RNA signatures as immunotherapy biomarkers using GSVA, ssGSEA, Weighted mean expression, and Specific algorithms (PredictIO).

// Compute signature scores
process GeneSigScore {
    tag "${study_id}"
    container 'bhklab/nextflow-env:latest'
    publishDir "${params.out_dir}/${study_id}", mode: 'copy'

    input:
    tuple val(study_id), path(signature_information), path(signature_data), path(icb_rda_path)

    output:
    tuple val(study_id), path("${study_id}_GeneSigScore.csv")

    script:
    """
    #!/usr/bin/env Rscript
    source('/R/load_libraries.R')

    signature <- read.csv("${signature_information}")
    signature\$Signature <- as.character(signature\$signature)
    signature\$method <- as.character(signature\$method)

    GeneSig_list <- list.files(path = '${signature_data}', pattern = '*.rda', full.names = TRUE)

    load("${icb_rda_path}")

    geneSig.score <- lapply(1:length(GeneSig_list), function(i) {
        load(GeneSig_list[i])
        sig_name <- substr(basename(GeneSig_list[i]), 1, nchar(basename(GeneSig_list[i])) - 4)

        method <- signature[signature\$Signature == sig_name, "method"]

        if (signature[signature\$Signature == sig_name, "method"] == "GSVA") {
            geneSig <- geneSigGSVA(dat.icb = ${study_id}, sig = sig, sig.name = sig_name, missing.perc = 0.5, const.int = 0.001, n.cutoff = 15, sig.perc = 0.8, study = "${study_id}")
            if (sum(!is.na(geneSig)) > 0) {
                geneSig <- geneSig[1,]
            }
        } else if (signature[signature\$Signature == sig_name, "method"] == "Weighted Mean") {
            geneSig <- geneSigMean(dat.icb = ${study_id}, sig = sig, sig.name = sig_name, missing.perc = 0.5, const.int = 0.001, n.cutoff = 15, sig.perc = 0.8, study = "${study_id}")
        } else if (signature[signature\$Signature == sig_name, "method"] == "ssGSEA") {
            geneSig <- geneSigssGSEA(dat.icb = ${study_id}, sig = sig, sig.name = sig_name, missing.perc = 0.5, const.int = 0.001, n.cutoff = 15, sig.perc = 0.8, study = "${study_id}")
            if (sum(!is.na(geneSig)) > 0) {
                geneSig <- geneSig[1,]
            }
        } else if (signature[signature\$Signature == sig_name, "method"] == "Specific Algorithm" & sig_name == "COX-IS_Bonavita") {
            geneSig <- geneSigCOX_IS(dat.icb = ${study_id}, sig = sig, sig.name = signature\$Signature[i], missing.perc = 0.5, const.int = 0.001, n.cutoff = 15, sig.perc = 0.8, study = "${study_id}")
        } else if (signature[signature\$Signature == sig_name, "method"] == "Specific Algorithm" & sig_name == "IPS_Charoentong") {
            geneSig <- geneSigIPS(dat.icb = ${study_id}, sig = sig, sig.name = signature\$Signature[i], missing.perc = 0.5, const.int = 0.001, n.cutoff = 15, sig.perc = 0.8, study = "${study_id}")
        } else if (signature[signature\$Signature == sig_name, "method"] == "Specific Algorithm" & sig_name == "PredictIO_Bareche") {
            geneSig <- geneSigPredictIO(dat.icb = ${study_id}, sig = sig, sig.name = signature\$Signature[i], missing.perc = 0.5, const.int = 0.001, n.cutoff = 15, sig.perc = 0.8, study = "${study_id}")
        } else if (signature[signature\$Signature == sig_name, "method"] == "Specific Algorithm" & sig_name == "IPRES_Hugo") {
            geneSig <- geneSigIPRES(dat.icb = ${study_id}, sig = sig, sig.name = signature\$Signature[i], missing.perc = 0.5, const.int = 0.001, n.cutoff = 15, sig.perc = 0.8, study = "${study_id}")
        } else if (signature[signature\$Signature == sig_name, "method"] == "Specific Algorithm" & sig_name == "PassON_Du") {
            geneSig <- geneSigPassON(dat.icb = ${study_id}, sig = sig, sig.name = signature\$Signature[i], missing.perc = 0.5, const.int = 0.001, n.cutoff = 15, sig.perc = 0.8, study = "${study_id}")
        } else if (signature[signature\$Signature == sig_name, "method"] == "Specific Algorithm" & sig_name == "IPSOV_Shen") {
            geneSig <- geneSigIPSOV(dat.icb = ${study_id}, sig = sig, sig.name = signature\$Signature[i], missing.perc = 0.5, const.int = 0.001, n.cutoff = 15, sig.perc = 0.8, study = "${study_id}")
        }

        if (sum(!is.na(geneSig)) > 0) {
            geneSig <- geneSig
        } else {
            geneSig <- rep(NA, ncol(${study_id}))
        }

        geneSig
    })

    geneSig.score <- do.call(rbind, geneSig.score)
    rownames(geneSig.score) <- substr(basename(GeneSig_list), 1, nchar(basename(GeneSig_list)) - 4)

    remove <- which(is.na(rowSums(geneSig.score)))
    if (length(remove) > 0) {
        geneSig.score <- geneSig.score[-remove, ]
    }
    write.csv(geneSig.score, file = "${study_id}_GeneSigScore.csv", row.names = TRUE)
    """
}


process GeneSig_AssociationOS {
    tag "${study_id}"
    container 'bhklab/nextflow-env:latest'
    publishDir "${params.out_dir}/${study_id}", mode: 'copy'

    input:
    tuple val(study_id), path(icb_rda_path), path(genescore_path), val(cancer_type), val(treatment)

    output:
    path("${study_id}_os_GeneSig_association.csv")


    script:
    """
    #!/usr/bin/env Rscript
    source('/R/load_libraries.R')

    load("${icb_rda_path}")
    geneSig.score <- read.csv("${genescore_path}", row.names = 1)

    res.all <- lapply(1:nrow(geneSig.score), function(k) {
        sig_name <- rownames(geneSig.score)[k]
        geneSig_vector <- as.numeric(geneSig.score[k, ])
        geneSig_vector <- geneSig_vector[!is.na(geneSig_vector)]

        res <- geneSigSurvCont(
            dat.icb = ${study_id},
            geneSig = geneSig_vector,
            time.censor = 36,
            n.cutoff = 15,
            study = "${study_id}",
            surv.outcome = "OS",
            sig.name = sig_name,
            cancer.type = "${cancer_type}",
            treatment = "${treatment}"
        )
        
        res
    })

    res.all <- do.call(rbind, res.all)
    res.all\$FDR <- p.adjust(res.all\$Pval, method="BH")
    res.all <- res.all[order(res.all\$FDR), ]

    write.csv(res.all, file = "${study_id}_os_GeneSig_association.csv", row.names = TRUE)
    """
}

process GeneSig_AssociationPFS {
    tag "${study_id}"
    container 'bhklab/nextflow-env:latest'
    publishDir "${params.out_dir}/${study_id}", mode: 'copy'

    input:
    tuple val(study_id), path(icb_rda_path), path(genescore_path), val(cancer_type), val(treatment)

    output:
    path("${study_id}_pfs_GeneSig_association.csv")

    script:
    """
    #!/usr/bin/env Rscript
    source('/R/load_libraries.R')

    load("${icb_rda_path}")
    
    geneSig.score <- read.csv("${genescore_path}", row.names = 1)

    res.all <- lapply(1:nrow(geneSig.score), function(k) {
        sig_name <- rownames(geneSig.score)[k]
        geneSig_vector <- as.numeric(geneSig.score[k, ])
        geneSig_vector <- geneSig_vector[!is.na(geneSig_vector)]
    
        res <- geneSigSurvCont(
            dat.icb = ${study_id},
            geneSig = geneSig_vector,  
            time.censor = 24,
            n.cutoff = 15,
            study = "${study_id}",
            surv.outcome = "PFS",
            sig.name = sig_name,
            cancer.type = "${cancer_type}",
            treatment = "${treatment}"
        )

        res
    })

    res.all <- do.call(rbind, res.all)
    res.all\$FDR <- p.adjust(res.all\$Pval, method="BH")
    res.all <- res.all[order(res.all\$FDR), ]

    write.csv(res.all, file = "${study_id}_pfs_GeneSig_association.csv", row.names = TRUE)
    """
}

// Repeat similar changes for other processes such as GeneSig_AssociationOS, GeneSig_AssociationResponse, etc.

process GeneSig_AssociationResponse {
    tag "${study_id}"
    container 'bhklab/nextflow-env:latest'
    publishDir "${params.out_dir}/${study_id}", mode: 'copy'

    input:
    tuple val(study_id), path(icb_rda_path), path(genescore_path), val(cancer_type), val(treatment)

    output:
    path("${study_id}_GeneSig_Response.csv")

    script:
    """
    #!/usr/bin/env Rscript
    source('/R/load_libraries.R')

    load("${icb_rda_path}")
    geneSig.score <- read.csv("${genescore_path}", row.names = 1)

    res.logreg <- lapply(1:nrow(geneSig.score), function(k){
        sig_name <- rownames(geneSig.score)[k]
        geneSig_vector <- as.numeric(geneSig.score[k, ])
        geneSig_vector <- geneSig_vector[!is.na(geneSig_vector)]

        res <- geneSigLogReg(dat.icb = ${study_id},
                            geneSig = geneSig_vector,
                            n.cutoff = 10,
                            study =  "${study_id}",
                            sig.name = sig_name,
                            n0.cutoff = 3, 
                            n1.cutoff = 3,
                            cancer.type = "${cancer_type}",
                            treatment = "${treatment}")
        
        res
    })

    res.logreg <- do.call(rbind, res.logreg)
    res.logreg\$FDR <- p.adjust(res.logreg\$Pval, method="BH")
    # Save as CSV file
    write.csv(res.logreg, file = "${study_id}_GeneSig_Response.csv", row.names = TRUE)
    """
}

/*
========================================================
SECTION: Gene level Meta Analysis 
========================================================
*/

/*
The following clinical multimodal immunotherapy datasets are publicly available on GitHub. These datasets are used in biomarker discovery for immunotherapy response through treatment-specific analyses.

Links:
- GitHub: https://github.com/bhklab/PredictioR/tree/main/data
*/

/*
-----------------------------------------------------------------------
SUBSECTION: Aggregating Associations through Meta-analysis (Pan-cancer)
-----------------------------------------------------------------------

Load public clinical multimodal immunotherapy datasets from GitHub or ORCESTRA 
for transparent biomarker discovery in immunotherapy response. For RNA profiles,
we use log2-transformed TPM data from protein-coding genes, filtering out genes 
with zero expression in at least 50% of samples. Only studies with at least 20 
patients are included.

Links:
- GitHub: https://github.com/bhklab/PredictioR/tree/main/data
- ORCESTRA: https://www.orcestra.ca/clinical_icb

Here './ICB_data' directory contains eight ICB data files for meta-analysis:
"ICB_Liu", "ICB_Padron", "ICB_Hugo", "ICB_Mariathasan", "ICB_Nathanson", 
"ICB_Riaz", "ICB_Miao", "ICB_Van_Allen"

Using gene CXCL9, to generalize the association with immunotherapy survival, 
we apply a meta-analysis approach to integrate findings across datasets for 
pan-cancer and per-cancer analysis.
*/
process MetaAnalysis_Gene_PanCancer {
    tag { "${params.gene_name} using ${io_outcome}" }
    container 'bhklab/nextflow-env:latest'
    publishDir "${params.out_dir}", mode: 'copy'

    input:
    tuple val(io_outcome), path(result_dir)
    // 'io_outcome' is immunotherapy outcome can be "OS", "PFS" or "Response"

    output:
    path "Meta_analysis_${io_outcome}_${params.gene_name}_pancancer.csv"

    script:
    """
    #!/usr/bin/env Rscript
    source('/R/load_libraries.R')

    # Determine the pattern based on io_outcome
    pattern <- if ('${io_outcome}' == 'OS') {
        '_cox_os.csv'
    } else if ('${io_outcome}' == 'PFS') {
        '_cox_pfs.csv'
    } else if ('${io_outcome}' == 'Response') {
        '_logregResponse.csv'
    } else {
        stop("Invalid io_outcome: ${io_outcome}. Must be 'OS', 'PFS', or 'Response'.")
    }

    # List all files that contain the pattern in their filenames within subdirectories
    result_files <- list.files(path = "${result_dir}", pattern = pattern, full.names = TRUE, recursive = TRUE)
    
    if (length(result_files) == 0) {
        stop("No files found matching the pattern. Please check the result directory and pattern.")
    }
    
    # Read each file and store the results
    res <- lapply(result_files, function(file) {
        df <- read.csv(file)
        df\$Study <- sub(pattern, "", basename(file))
        
        
        df
    })

    # Combine the results into a single data frame
    assoc.res <- do.call(rbind, res)
    assoc.res <- assoc.res[!is.na(assoc.res\$Coef), ]

    # Adjust p-values for multiple comparisons using the Benjamini-Hochberg method
    assoc.res\$FDR <- p.adjust(assoc.res\$Pval, method = "BH")

    # Meta-analysis for a gene across datasets
    res_meta_pancancer <- metafun(
        coef = assoc.res\$Coef, 
        se = assoc.res\$SE, 
        study = assoc.res\$Study, 
        pval = assoc.res\$Pval, 
        n = assoc.res\$N, 
        cancer.type = assoc.res\$Cancer_type, 
        treatment = assoc.res\$Treatment, 
        feature = "${params.gene_name}", 
        cancer.spec = FALSE, 
        treatment.spec = FALSE
    )

    # Save the results to a CSV file
    write.csv(data.frame(res_meta_pancancer), file = "Meta_analysis_${io_outcome}_${params.gene_name}_pancancer.csv", row.names = FALSE)
    """
}

/*
-----------------------------------------------------------------------
SUBSECTION: Aggregating Associations through Meta-analysis (Per-cancer)
-----------------------------------------------------------------------

For cancer-specific analysis, consider meta-analysis when there are at least 3 datasets.
*/

process MetaAnalysis_Gene_PerCancer {
    tag { "${params.gene_name} using ${io_outcome}" }
    container 'bhklab/nextflow-env:latest'
    publishDir "${params.out_dir}", mode: 'copy'

    input:
    tuple val(io_outcome), path(result_dir)
    // 'io_outcome' is immunotherapy outcome can be "OS", "PFS" or "Response"

    output:
    path "Meta_analysis_${io_outcome}_${params.gene_name}_percancer.csv"

    script:
    """
    #!/usr/bin/env Rscript
    source('/R/load_libraries.R')

    # Determine the pattern based on io_outcome
    pattern <- if ('${io_outcome}' == 'OS') {
        '_cox_os.csv'
    } else if ('${io_outcome}' == 'PFS') {
        '_cox_pfs.csv'
    } else if ('${io_outcome}' == 'Response') {
        '_logregResponse.csv'
    } else {
        stop("Invalid io_outcome: ${io_outcome}. Must be 'OS', 'PFS', or 'Response'.")
    }

    # List all files that contain the pattern in their filenames within subdirectories
    result_files <- list.files(path = "${result_dir}", pattern = pattern, full.names = TRUE, recursive = TRUE)
    
    
    if (length(result_files) == 0) {
        stop("No files found matching the pattern. Please check the result directory and pattern.")
    }
    
    # Read each file and store the results
    res <- lapply(result_files, function(file) {
        df <- read.csv(file)
        df\$Study <- sub(pattern, "", basename(file))
        
        df
    })

    # Combine the results into a single data frame
    assoc.res <- do.call(rbind, res)
    assoc.res <- assoc.res[!is.na(assoc.res\$Coef), ]

    # Check for empty assoc.res
    if (nrow(assoc.res) == 0) {
        stop("No valid rows in assoc.res data frame. Please check the input files.")
    }

    # Adjust p-values for multiple comparisons using the Benjamini-Hochberg method
    assoc.res\$FDR <- p.adjust(assoc.res\$Pval, method = "BH")

    # Meta-analysis for a gene across datasets
    res_meta_pancancer <- metafun(
        coef = assoc.res\$Coef, 
        se = assoc.res\$SE, 
        study = assoc.res\$Study, 
        pval = assoc.res\$Pval, 
        n = assoc.res\$N, 
        cancer.type = assoc.res\$Cancer_type, 
        treatment = assoc.res\$Treatment, 
        feature = "${params.gene_name}", 
        cancer.spec = FALSE, 
        treatment.spec = FALSE
    )

    # Treatment-specific meta-analysis for a gene across datasets
    res_meta_percancer <- metaPerCanfun(coef = assoc.res\$Coef, se = assoc.res\$SE, study = assoc.res\$Study, pval = assoc.res\$Pval, n = assoc.res\$N, cancer.type = assoc.res\$Cancer_type, treatment = assoc.res\$Treatment, feature = "${params.gene_name}", cancer.spec = TRUE)

    # Combine all meta_summery results into a single data frame
    meta_summery_combined <- do.call(rbind, lapply(res_meta_percancer, function(x) x\$meta_summery))
    write.csv(meta_summery_combined, file = "Meta_analysis_${io_outcome}_${params.gene_name}_percancer.csv", row.names = FALSE)
    """
}

/*
--------------------------------------------------------
Sigevel Meta-analysis : Pan-cancer 
--------------------------------------------------------
*/

/*
--------------------------------------------------------
Sigevel Meta-analysis : Pan-cancer 
--------------------------------------------------------
*/

process MetaAnalysis_Sig_PanCancer {
    tag " sigGenes using ${io_outcome}"
    container 'bhklab/nextflow-env:latest'
    publishDir "${params.out_dir}", mode: 'copy'

    input:
    tuple val(io_outcome), path(result_dir) 
    
    output:
    path "Meta_analysis_Sig_${io_outcome}_pancancer.csv"

    script:
    """
    #!/usr/bin/env Rscript
    source('/R/load_libraries.R')

    # Define the directory where your files are located
    sig_level_result_dir <- '${result_dir}'

    # Determine the pattern based on io_outcome
    pattern <- ifelse('${io_outcome}' == 'OS', '_os_GeneSig', ifelse('${io_outcome}' == 'PFS', '_pfs_GeneSig', 'Response_GeneSig'))

    # List all files that contain the pattern in their filenames
    sig_files <- list.files(path = sig_level_result_dir, pattern = pattern, full.names = TRUE, recursive = TRUE)

    # Read each file and store the results
    res <- lapply(sig_files, function(file) {
        df <- read.csv(file)
        df
    })

    # Combine the results into a single data frame
    res <- do.call(rbind, res)
    res <- res[!is.na(res\$Coef), ]

    # Convert to data frame
    df <- res
    signature <- unique(df\$Gene)

    # Perform meta-analysis on each gene signature
    AllGeneSig_meta <- lapply(1:length(signature), function(j) {
        res <- metafun(coef = df[df\$Gene == signature[j], "Coef"],
                       se = df[df\$Gene == signature[j], "SE"],
                       study  = df[df\$Gene == signature[j], "Study"],
                       pval = df[df\$Gene == signature[j], "Pval"],
                       n = df[df\$Gene == signature[j], "N"],
                       cancer.type = df[df\$Gene == signature[j], "Cancer_type"],
                       treatment = df[df\$Gene == signature[j], "Treatment"],
                       cancer.spec = FALSE,
                       treatment.spec = FALSE,
                       feature = unique(df[df\$Gene == signature[j], "Gene"]))

        res\$meta_summery
    })

    # Combine meta-analysis results
    AllGeneSig_meta <- do.call(rbind, AllGeneSig_meta)
    AllGeneSig_meta <- AllGeneSig_meta[!is.na(AllGeneSig_meta\$Coef), ]
    AllGeneSig_meta\$FDR <- p.adjust(AllGeneSig_meta\$Pval, method = "BH")
    AllGeneSig_meta <- AllGeneSig_meta[order(AllGeneSig_meta\$FDR), ]

    # Save the results to a CSV file
    write.csv(AllGeneSig_meta, file = "Meta_analysis_Sig_${io_outcome}_pancancer.csv", row.names = FALSE)
    """
}


/*
--------------------------------------------------------
Sig-level Meta-analysis : Per-cancer 
--------------------------------------------------------
*/

process MetaAnalysis_Sig_PerCancer {
    tag " sigGenes using ${io_outcome}"
    container 'bhklab/nextflow-env:latest'
    publishDir "${params.out_dir}", mode: 'copy'

    input:
    tuple val(io_outcome), path(result_dir) 
    
    output:
    path "Meta_analysis_Sig_${io_outcome}_percancer.csv"

    script:
    """
    #!/usr/bin/env Rscript
    source('/R/load_libraries.R')

    # Define the directory where your files are located
    sig_level_result_dir <- '${result_dir}'

    # Determine the pattern based on io_outcome
    pattern <- ifelse('${io_outcome}' == 'OS', '_os_GeneSig', ifelse('${io_outcome}' == 'PFS', '_pfs_GeneSig', 'GeneSig_Response'))

    # List all files that contain the pattern in their filenames
    sig_files <- list.files(path = sig_level_result_dir, pattern = pattern, full.names = TRUE, recursive = TRUE)

    # Read each file and store the results
    res <- lapply(sig_files, function(file) {
        df <- read.csv(file)
        df
    })

    # Combine the results into a single data frame
    res <- do.call(rbind, res)
    res <- res[!is.na(res\$Coef), ]

    # Convert to data frame
    df <- res
    signature <- unique(df\$Gene)

    # Perform per-cancer meta-analysis on each gene signature
    AllGeneSig_meta <- lapply(1:length(signature), function(j) {
        sub_df <- df[df\$Gene == signature[j], ]
        if (nrow(sub_df) >= 3) {
            res <- metaPerCanfun(coef = sub_df\$Coef,
                                 se = sub_df\$SE,
                                 study  = sub_df\$Study,
                                 pval = sub_df\$Pval,
                                 n = sub_df\$N,
                                 cancer.type = sub_df\$Cancer_type,
                                 treatment = sub_df\$Treatment,
                                 cancer.spec = TRUE,
                                 feature = unique(sub_df\$Gene))

            percan_res <- lapply(1:length(res), function(i) {
                res[[i]]\$meta_summery
            })

            percan_res <- do.call(rbind, percan_res)
        } else {
            percan_res <- data.frame(Cancer_type = "Not Applicable",
                                     Gene = signature[j],
                                     Coef = NA,
                                     SE = NA,
                                     CI_lower = NA,
                                     CI_upper = NA,
                                     Pval = NA,
                                     I2 = NA,
                                     Q_Pval = NA)
        }
        percan_res
    })

    AllGeneSig_meta <- do.call(rbind, AllGeneSig_meta)
    AllGeneSig_meta <- AllGeneSig_meta[!is.na(AllGeneSig_meta\$Coef), ]

    # FDR adjustment
    group <- unique(AllGeneSig_meta\$Cancer_type)
    AllGeneSig_meta <- lapply(1:length(group), function(k) {
        sub_df <- AllGeneSig_meta[AllGeneSig_meta\$Cancer_type == group[k], ]
        sub_df\$FDR <- p.adjust(sub_df\$Pval, method = "BH")
        sub_df
    })

    AllGeneSig_meta <- do.call(rbind, AllGeneSig_meta)
    AllGeneSig_meta <- AllGeneSig_meta[order(AllGeneSig_meta\$FDR), ]

    # Save the results to a CSV file
    write.csv(AllGeneSig_meta, file = "Meta_analysis_Sig_${io_outcome}_percancer.csv", row.names = FALSE)
    """
}

workflow {

    /*
    ===============================================================
    SECTION: Load RDA Data and Extract Expression and Clinical Data
    ===============================================================
    */

    // List all .rda files in the data directory
    icb_rda_files = Channel.fromPath("${params.icb_data_dir}/*.rda")

    // Load and extract data from each .rda file
    extracted_data = icb_rda_files.map { file ->
        study_id = file.baseName
        tuple(study_id, file)
    } | LoadAndExtractData

    // Map the extracted data with cancer type, treatment information, and genes
    extracted_data_with_info = extracted_data.map { study_id, expr_file, clin_file, annot_file ->
        cancer_type = cancer_type_map[study_id]
        treatment = treatment_map[study_id]

        // Choose genes as a single gene or a vector of genes
        genes = 'c("CXCL9", "CXCL10", "TIGIT", "CD83", "STAT1", "CXCL11", "CXCL13", "CD8A", "CTLA4")'
        tuple(study_id, expr_file, clin_file, cancer_type, treatment, genes)
    }
    

    /*
    ========================================================
    SECTION: Gene Association Analysis
    ========================================================
    */

    // OS analysis
    geneassosiation_os_results = extracted_data_with_info | GeneAssociationOS

    // PFS analysis
    geneassosiation_pfs_results = extracted_data_with_info | GeneAssociationPFS

    // Immunotherapy response analysis (R vs NR)
    geneassosiation_response_results = extracted_data_with_info | GeneAssociationResponse

 
    /*
    ========================================================
    SECTION: Signature Score Computation 
    ========================================================
    */

    // Example process using the icb_rda_files channel
    icb_rda_files.map { file -> tuple(file.baseName, file) }.set { query_ch }

    // Signature information and data
    sigs_info_path = file("${params.sig_summery_dir}/signature_information.csv")
    signature_data = file(params.sig_data_dir)

    signature_analysis = icb_rda_files.map { file ->
        def study_id = file.baseName
        tuple(study_id, sigs_info_path, signature_data, file)
    } | GeneSigScore


    signature_analysis_with_info = signature_analysis.map { study_id, genescore_path ->
        cancer_type = cancer_type_map[study_id]
        treatment = treatment_map[study_id]
        rda_path = file("${params.icb_data_dir}/${study_id}.rda")
        genescore_full_path = file("${params.out_dir}/${study_id}/${study_id}_GeneSigScore.csv")
        tuple(study_id, rda_path, genescore_full_path, cancer_type, treatment)
    }
    
    /*
    ========================================================
    SECTION: Signature Level Analysis
    ========================================================
    */

    // OS analysis for signatures
    signature_os_results = signature_analysis_with_info | GeneSig_AssociationOS

    // PFS analysis for signatures
    signature_pfs_results = signature_analysis_with_info | GeneSig_AssociationPFS

    // Response analysis for signatures
    signature_response_results = signature_analysis_with_info | GeneSig_AssociationResponse

/* 
After you use `nextflow run main.nf`, uncomment this section and use `nextflow run main.nf -resume` to have the meta-analysis section added
*/

/*
========================================================
Meta Analysis for Genelevel + SigLevel
=======================================================

// 1. Pan-cancer
// "OS", "PFS" or "Response" can be used as io_outcome
// you can choose "OS", "PFS" or "Response" as io_outcome
gene_os_result_files = Channel.of(['OS', file("${params.out_dir}")])
gene_os_result_files | MetaAnalysis_Gene_PanCancer

// 2. Per-cancer
gene_response_result_files = Channel.of(['Response', file("${params.out_dir}")])
gene_response_result_files | MetaAnalysis_Gene_PerCancer

// 1. Pan-cancer
// you can choose "OS", "PFS" or "Response" as io_outcome
meta_analysis_sig_pancancer = Channel.of(['OS', file("${params.out_dir}")])
meta_analysis_sig_pancancer | MetaAnalysis_Sig_PanCancer

// 2. Per-cancer
meta_analysis_sig_percancer = Channel.of(['Response', file("${params.out_dir}")])
meta_analysis_sig_percancer | MetaAnalysis_Sig_PerCancer
*/

}
