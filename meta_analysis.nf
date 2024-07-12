#!/usr/bin/env nextflow
nextflow.enable.dsl=2

params.icb_data_dir = './ICB_data'
params.out_dir = './output/meta_analysis_output'
params.sig_level_result_dir = './output/signatures_level_output'


// Set gene of interest
params.gene_name = "CXCL9"

// Use all .rda files from './ICB_data'.
// Set the cancer type and treatment for each dataset respectively.
params.cancer_types = '["Melanoma", "Pancreas", "Melanoma", "Bladder", "Melanoma", "Melanoma", "Kidney", "Melanoma"]'
params.treatment_types = '["PD-1/PD-L1", "PD-1/PD-L1", "PD-1/PD-L1", "PD-1/PD-L1", "CTLA4", "IO+combo", "PD-1/PD-L1", "CTLA4"]'

log.info """
P R E D I C T I O - N F   P I P E L I N E (for Gene Level + Signature Level)
========================================================================
ICB Data Directory: ${params.icb_data_dir}
Output Directory  : ${params.out_dir}
Cancer Types      : ${params.cancer_types}
Treatments        : ${params.treatment_types}
Genenames         : ${params.gene_name}
""".stripIndent()

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
    container 'nextflow-env:latest'
    publishDir "${params.out_dir}", mode: 'copy'

    input:
    path rda_files_path
    val io_outcome // immunotherapy outcome can be "OS", "PFS" or "Response"

    output:
    path "Meta_analysis_${io_outcome}_${params.gene_name}_pancancer.csv"

    script:
    """
    #!/usr/bin/env Rscript
    source('/R/load_libraries.R')

    # Load all .rda files 
    lapply(list.files(path = '${rda_files_path}', pattern = '*.rda', full.names = TRUE), function(file) { load(file, envir = .GlobalEnv) })

    # Create a list of the loaded objects using their actual names
    loaded_objects <- basename(list.files(path = '${rda_files_path}', pattern = '*.rda', full.names = TRUE))
    loaded_objects <- substr(loaded_objects, 1, nchar(loaded_objects) - 4) # remove .rda from name
    expr <- mget(loaded_objects, envir = .GlobalEnv)

    # Define the cancer types and treatment types vectors
    cancer_types <- fromJSON('${params.cancer_types}')

    treatment_types <- fromJSON('${params.treatment_types}')

    # Apply a function over the loaded datasets to perform survival or response analysis
    assoc.res <- lapply(1:length(expr), function(k){
        if ('${io_outcome}' == 'OS') {
            geneSurvCont(dat.icb = expr[[k]], time.censor = 36, missing.perc = 0.5, const.int = 0.001, n.cutoff = 15, feature = "${params.gene_name}", study = names(expr)[k], surv.outcome = '${io_outcome}', cancer.type = cancer_types[k], treatment = treatment_types[k])
        } else if ('${io_outcome}' == 'PFS') {
            geneSurvCont(dat.icb = expr[[k]], time.censor = 24, missing.perc = 0.5, const.int = 0.001, n.cutoff = 15, feature = "${params.gene_name}", study = names(expr)[k], surv.outcome = '${io_outcome}', cancer.type = cancer_types[k], treatment = treatment_types[k])
        } else if ('${io_outcome}' == 'Response') {
            geneLogReg(dat.icb = expr[[k]], missing.perc = 0.5, const.int = 0.001, n.cutoff = 15, feature = "${params.gene_name}", study = names(expr)[k], n0.cutoff = 3, n1.cutoff = 3, cancer.type = cancer_types[k], treatment = treatment_types[k])
        } else {
            print(paste("You can only use 'OS', 'PFS', or 'Response' as '${io_outcome}' input."))
            return(NULL)
        }
    })

    # Combine results into one data frame
    assoc.res <- do.call(rbind, assoc.res)

    # Adjust p-values for multiple comparisons using the Benjamini-Hochberg method
    assoc.res\$FDR <- p.adjust(assoc.res\$Pval, method = "BH")

    # Meta-analysis for a gene across datasets
    res_meta_pancancer <- metafun(coef = assoc.res\$Coef, se = assoc.res\$SE, study = assoc.res\$Study, pval = assoc.res\$Pval, n = assoc.res\$N, cancer.type = assoc.res\$Cancer_type, treatment = assoc.res\$Treatment, feature = "${params.gene_name}", cancer.spec = FALSE, treatment.spec = FALSE)

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
    container 'nextflow-env:latest'
    publishDir "${params.out_dir}", mode: 'copy'

    input:
    path rda_files_path
    val io_outcome // immunotherapy outcome can be "OS", "PFS" or "Response"

    output:
    path "Meta_analysis_${io_outcome}_${params.gene_name}_percancer.csv"

    script:
    """
    #!/usr/bin/env Rscript
    source('/R/load_libraries.R')

    # Load all .rda files 
    lapply(list.files(path = '${rda_files_path}', pattern = '*.rda', full.names = TRUE), function(file) { 
        load(file, envir = .GlobalEnv)
    })

    # Create a list of the loaded objects using their actual names
    loaded_objects <- basename(list.files(path = '${rda_files_path}', pattern = '*.rda', full.names = TRUE))
    loaded_objects <- substr(loaded_objects, 1, nchar(loaded_objects) - 4) # remove .rda from name
    expr <- mget(loaded_objects, envir = .GlobalEnv)

    # Define the cancer types and treatment types vectors
    cancer_types <- fromJSON('${params.cancer_types}')
    treatment_types <- fromJSON('${params.treatment_types}')

    # Apply a function over the loaded datasets to perform survival or response analysis
    assoc.res <- lapply(1:length(expr), function(k){
        print(paste("Processing dataset:", names(expr)[k])) # Debug
        if ('${io_outcome}' == 'OS') {
            return(geneSurvCont(dat.icb = expr[[k]], time.censor = 36, missing.perc = 0.5, const.int = 0.001, n.cutoff = 15, feature = "${params.gene_name}", study = names(expr)[k], surv.outcome = '${io_outcome}', cancer.type = cancer_types[k], treatment = treatment_types[k]))
        } else if ('${io_outcome}' == 'PFS') {
            return(geneSurvCont(dat.icb = expr[[k]], time.censor = 24, missing.perc = 0.5, const.int = 0.001, n.cutoff = 15, feature = "${params.gene_name}", study = names(expr)[k], surv.outcome = '${io_outcome}', cancer.type = cancer_types[k], treatment = treatment_types[k]))
        } else if ('${io_outcome}' == 'Response') {
            return(geneLogReg(dat.icb = expr[[k]], missing.perc = 0.5, const.int = 0.001, n.cutoff = 15, feature = "${params.gene_name}", study = names(expr)[k], n0.cutoff = 3, n1.cutoff = 3, cancer.type = cancer_types[k], treatment = treatment_types[k]))
        } else {
            print(paste("You can only use 'OS', 'PFS', or 'Response' as '${io_outcome}' input."))
            return(NULL)
        }
    })

    # Combine results into one data frame
    assoc.res <- do.call(rbind, assoc.res)

    # Adjust p-values for multiple comparisons using the Benjamini-Hochberg method
    assoc.res\$FDR <- p.adjust(assoc.res\$Pval, method = "BH")

    # Treatment-specific meta-analysis for a gene across datasets
    res_meta_percancer <- metaPerCanfun(coef = assoc.res\$Coef, se = assoc.res\$SE, study = assoc.res\$Study, pval = assoc.res\$Pval, n = assoc.res\$N, cancer.type = assoc.res\$Cancer_type, treatment = assoc.res\$Treatment, feature = "${params.gene_name}", cancer.spec = TRUE)

    # Combine all meta_summery results into a single data frame
    meta_summery_combined <- do.call(rbind, lapply(res_meta_percancer, function(x) x\$meta_summery))
    write.csv(meta_summery_combined, file = "Meta_analysis_${io_outcome}_${params.gene_name}_percancer.csv", row.names = FALSE)
    """
}

/*
========================================================
SECTION: Signature level Meta Analysis 
========================================================
*/  

/*
----------------------------
SUBSECTION: Pan-cancer
----------------------------
Uses the _os, _pfs, and _Response files from signature_level_analysis.nf using
all ICB datasets in ./ICB_data. Results are stored in sig_level_result_dir.

Note: You should have at least three results (studies) to proceed.
*/


process MetaAnalysis_Sig_PanCancer {
    tag " sigGenes using ${io_outcome}"
    container 'nextflow-env:latest'
    publishDir "${params.out_dir}", mode: 'copy'

    input:
    path sig_level_result_dir 
    val io_outcome // can be 'OS', 'PFS', or 'Response' (RvsNR)

    output:
    path "Meta_analysis_Sig_${io_outcome}_pancancer.csv"

    script:
    """
    #!/usr/bin/env Rscript
    source('/R/load_libraries.R')

    # Define the directory where your files are located
    sig_level_result_dir <- '${sig_level_result_dir}'

    # Determine the pattern based on io_outcome
    pattern <- ifelse('${io_outcome}' == 'OS', '_os', ifelse('${io_outcome}' == 'PFS', '_pfs', '_Response'))

    # List all files that contain the pattern in their filenames
    sig_files <- list.files(path = sig_level_result_dir, pattern = pattern, full.names = TRUE)

    # Read each file and store the results
    res <- lapply(sig_files, function(file) {
        read.csv(file)
    })

    # Combine the results into a single data frame
    res <- do.call(rbind, res)
    res <- res[!is.na(res\$Coef), ]

    # Convert to data frame
    df <- res
    signature <- unique(df\$Gene)

    # Perform meta-analysis on each gene signature
    AllGeneSig_meta <- lapply(1:length(signature), function(j) {
        print(j)
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
----------------------------
SUBSECTION: Per-cancer)
----------------------------
It uses the _os, _pfs, and _Response files from signature_level_analysis.nf using
all ICB datasets in ./ICB_data. Results are stored in sig_level_result_dir.

*/

process MetaAnalysis_Sig_PerCancer {
    tag "sigGenes using ${io_outcome}"
    container 'nextflow-env:latest'
    publishDir "${params.out_dir}", mode: 'copy'

    input:
    path sig_level_result_dir
    val io_outcome 

    output:
    path "Meta_analysis_Sig_${io_outcome}_percancer.csv"

    script:
    """
    #!/usr/bin/env Rscript
    source('/R/load_libraries.R')

    # Define the directory where your files are located
    sig_level_result_dir <- '${sig_level_result_dir}'

    # Determine the pattern based on io_outcome
    pattern <- ifelse('${io_outcome}' == 'OS', '_os', ifelse('${io_outcome}' == 'PFS', '_pfs', '_Response'))

    # List all files that contain the pattern in their filenames
    sig_files <- list.files(path = sig_level_result_dir, pattern = pattern, full.names = TRUE)

    # Read each file and store the results
    res <- lapply(sig_files, function(file) {
        read.csv(file)
    })

    # Combine the results into a single data frame
    res <- do.call(rbind, res)
    res <- res[!is.na(res\$Coef), ]

    # Convert to data frame
    df <- res
    signature <- unique(df\$Gene)

    # Perform per-cancer meta-analysis on each gene signature
    AllGeneSig_meta <- lapply(1:length(signature), function(j) {
        print(j)
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
    =================================
    SECTION: Gene level Meta Analysis 
    =================================
    */

    // Gene-level Meta Analysis
    all_icb_data = file(params.icb_data_dir)

    /*
    Pan-Cancer
    */

    // Perform meta-analysis for OS; you can use outcome "PFS" or "Response"
    MetaAnalysis_Gene_PanCancer(all_icb_data, io_outcome = "OS")

    // Uncomment below lines to handle each outcome type individually
    //MetaAnalysis_Gene_PanCancer(all_icb_data, io_outcome = "PFS")
    //MetaAnalysis_Gene_PanCancer(all_icb_data, io_outcome = "Response")

    /*
    Per-Cancer
    */

    MetaAnalysis_Gene_PerCancer(all_icb_data, io_outcome = "Response")

    // Uncomment below lines to handle each outcome type individually
    // MetaAnalysis_Gene_PerCancer(all_icb_data, io_outcome = "OS")
    //MetaAnalysis_Gene_PerCancer(all_icb_data, io_outcome = "PFS")

    /*
    =================================
    SECTION: Signature level Meta Analysis 
    =================================
    */

    /*
    Pan-Cancer
    */

    sig_level_results = file(params.sig_level_result_dir)

    //MetaAnalysis_Sig_PanCancer(sig_level_results, io_outcome = "OS" )

    // Uncomment below lines to handle each outcome type individually
    //MetaAnalysis_Sig_PanCancer(sig_level_results, io_outcome = "PFS")
    //MetaAnalysis_Sig_PanCancer(sig_level_results, io_outcome = "Response")

    /*
    Per-Cancer
    */

    sig_level_results = file(params.sig_level_result_dir)
    MetaAnalysis_Sig_PerCancer(sig_level_results, io_outcome = "OS" )

    // Uncomment below lines to handle each outcome type individually
    //MetaAnalysis_Sig_PerCancer(sig_level_results, io_outcome = "PFS")
    //MetaAnalysis_Sig_PerCancer(sig_level_results, io_outcome = "Response")

    }