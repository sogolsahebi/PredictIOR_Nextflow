#!/usr/bin/env nextflow
nextflow.enable.dsl=2

params.study_id = 'ICB_small_Mariathasan' // Set study_id (name of the .rda file in ./ICB_data directory)
params.data_dir = './ICB_data'
params.out_dir = './output'
params.cancer_type = 'Bladder'
params.treatment = 'PD-1/PD-L1'

params.cancer_types = '["Melanoma", "Pancreas", "Melanoma", "Bladder", "Melanoma", "Melanoma", "Kidney", "Melanoma"]'
params.treatment_types = '["PD-1/PD-L1", "PD-1/PD-L1", "PD-1/PD-L1", "PD-1/PD-L1", "CTLA4", "IO+combo", "PD-1/PD-L1", "CTLA4"]'
params.gene_name = "CXCL9"

log.info """
P R E D I C T I O - N F   P I P E L I N E (Gene Level Analsyis )
================================================================
Study ID             : ${params.study_id}
ICB Data Directory   : ${params.data_dir}
Output Directory     : ${params.out_dir}
Cancer Type          : ${params.cancer_type}
Treatment            : ${params.treatment}
""".stripIndent()

/*
========================================================
SECTION: Load Immunotherapy Datasets
========================================================
*/

/*
Load public clinical multimodal immunotherapy datasets from GitHub or ORCESTRA for transparent biomarker discovery in immunotherapy response. For RNA profiles, we use log2-transformed TPM data from protein-coding genes, filtering out genes with zero expression in at least 50% of samples. Only studies with at least 20 patients are included.

For example, the Mariathasan dataset (PMID 29443960) includes RNA expression, clinical data, and gene metadata for 195 patients with 17,993 protein-coding genes, focused on bladder cancer and PD-1/PD-L1 treatment.

Links:
- GitHub: https://github.com/bhklab/PredictioR/tree/main/data
- ORCESTRA: https://www.orcestra.ca/clinical_icb

The detailed clinical characteristics of the Mariathasan dataset include cancer type, age, sex, response (R vs. NR), and overall survival (OS).
*/

// Process to load RDA data and extract expression, clinical data, and annotation data
process LoadAndExtractData {
    tag "${params.study_id}"
    container 'nextflow-env:latest'
    publishDir "${params.out_dir}", mode: 'copy'

    input:
    path rda_file

    output:
    path "${params.study_id}_expr.csv"
    path "${params.study_id}_clin.csv"
    path "${params.study_id}_annot.csv"

    script:
    """
    #!/usr/bin/env Rscript
    source('/R/load_libraries.R')

    load("${rda_file}")

    expr <- assay(${params.study_id})
    clin <- as.data.frame(colData(${params.study_id}))
    annot <- as.data.frame(rowData(${params.study_id}))
    
    write.csv(expr, "${params.study_id}_expr.csv", row.names = TRUE)
    write.csv(clin, "${params.study_id}_clin.csv", row.names = FALSE)
    write.csv(annot, "${params.study_id}_annot.csv", row.names = TRUE)
    """
}

/*
========================================================
SECTION: Biomarkers and Immunotherapy Response Association - Gene Association Analysis

These section functions you can either:

Option 1: Use dat.icb = expr (data frame) with clin = clin (data frame) for clinical data
Option 2: Load the RDA file with load("${rda_file}"), then use  dat.icb = '${params.study_id}' (SummarizedExperiment object) with clin = NULL
========================================================
*/

/*
Assessing the association of specific biomarkers with:
- Immunotherapy response (R vs NR)
- Immunotherapy survival (for time-to-event analysis: OS and PFS)

P-values are corrected for multiple testing using the Benjamini-Hochberg (FDR) method, with significance set at p-values or FDR ≤ 5%.

Note: Coef represents the log hazard ratio (logHR) or log odds ratio (logOR) depending on the analysis.
*/

/*
--------------------------------------------------------------------
SUBSECTION: Overall survival (OS) or Progression-free survival (PFS)
--------------------------------------------------------------------
*/


process GeneAssociationSurvival {
    tag "${params.study_id}"
    container 'nextflow-env:latest'
    publishDir "${params.out_dir}", mode: 'copy'

    input:
    path expr_file
    path clin_file
    val survival_type
    val genes

    output:
    path("${params.study_id}_cox_${survival_type}_genes.csv")

    script:
    """
    #!/usr/bin/env Rscript
    source('/R/load_libraries.R')
    
    expr <- read.csv("${expr_file}", row.names = 1)
    clin <- read.csv("${clin_file}")

    cox_result <- geneSurvCont(
        dat.icb = expr,
        clin = clin,
        time.censor = 36,
        missing.perc = 0.5,
        const.int = 0.001,
        n.cutoff = 15,
        feature = ${genes},
        study = "${params.study_id}",
        surv.outcome = "${survival_type}",
        cancer.type = "${params.cancer_type}",
        treatment = "${params.treatment}"
    )

    # Adjust p-values for multiple testing
    cox_result\$FDR <- p.adjust(cox_result\$Pval, method = "BH")
    cox_result <- cox_result[order(cox_result\$FDR), ]
    write.csv(cox_result, file = "${params.study_id}_cox_${survival_type}_genes.csv", row.names = FALSE)
    """
}

/*
--------------------------------------------------------
SUBSECTION: Immunotherapy response (R vs NR)
--------------------------------------------------------
*/

process GeneAssociationResponse {
    tag "${params.study_id}"
    container 'nextflow-env:latest'
    publishDir "${params.out_dir}", mode: 'copy'

    input:
    path expr_file
    path clin_file
    val genes

    output:
    path("${params.study_id}_Response.csv")

    script:
    """
    #!/usr/bin/env Rscript
    source('/R/load_libraries.R')

    expr <- read.csv("${expr_file}", row.names = 1)
    clin <- read.csv("${clin_file}")

    logreg <- geneLogReg(
        dat.icb = expr,
        clin = clin,
        missing.perc = 0.5,
        const.int = 0.001,
        n.cutoff = 15,
        feature = ${genes},
        study = "${params.study_id}", 
        n0.cutoff = 5, 
        n1.cutoff = 5,
        cancer.type = "${params.cancer_type}",
        treatment = "${params.treatment}"
    )

    # Adjust P-values and sort by FDR
    logreg <- logreg[order(logreg\$FDR <- p.adjust(logreg\$Pval, method = "BH"))]

    # Save as CSV file
    write.csv(logreg, file = "${params.study_id}_Response.csv", row.names = FALSE)
    """
}

/*
========================================================
SECTION: Meta Analysis Section
========================================================
*/

/*
The following clinical multimodal immunotherapy datasets are publicly available on the Github. These datasets are used in biomarker discovery for immunotherapy response through treatment-specific analyses

Links:
- GitHub: https://github.com/bhklab/PredictioR/tree/main/data
*/

// Set the cancer type and treatment for each dataset respectively.
params.cancer_types = '["Melanoma", "Pancreas", "Melanoma", "Bladder", "Melanoma", "Melanoma", "Kidney", "Melanoma"]'
params.treatment_types = '["PD-1/PD-L1", "PD-1/PD-L1", "PD-1/PD-L1", "PD-1/PD-L1", "CTLA4", "IO+combo", "PD-1/PD-L1", "CTLA4"]'
params.gene_name = "CXCL9"


/*
-----------------------------------------------------------------------
SUBSECTION: Aggregating Associations through Meta-analysis (Per-treatment)
-----------------------------------------------------------------------
*/

process MetaAnalysis_pertreatment {
    tag {params.gene_name}
    container 'nextflow-env:latest'
    publishDir "${params.out_dir}", mode: 'copy'

    input:
    path rda_files_path
    val survival_type

    output:
    path "Meta_analysis_${survival_type}_${params.gene_name}_pertreatment.csv"

    script:
    """
    #!/usr/bin/env Rscript
    source('/R/load_libraries.R')

    # Load all .rda files into the global environment
    lapply(list.files(path = '${rda_files_path}', pattern = '*.rda', full.names = TRUE), function(file) {
        load(file, envir = .GlobalEnv)
    })

    # Create a list of the loaded objects using their actual names
    loaded_objects <- ls(pattern = "^ICB_")
    expr <- mget(loaded_objects, envir = .GlobalEnv)

    # Define the cancer types and treatment types vectors
    cancer_types <- fromJSON('${params.cancer_types}')
    treatment_types <- fromJSON('${params.treatment_types}')

    if (length(expr) != length(cancer_types) || length(expr) != length(treatment_types)) {
        stop("Mismatch in the length of loaded objects, cancer types, and treatment types.")
    }

    # Apply a function over the loaded datasets to perform survival or response analysis
    assoc.res <- lapply(1:length(expr), function(k){
        if ('${survival_type}' == 'OS' || '${survival_type}' == 'PFS') {
            geneSurvCont(
                dat.icb = expr[[k]],
                time.censor = 36,
                missing.perc = 0.5,
                const.int = 0.001,
                n.cutoff = 15,
                feature = "${params.gene_name}",
                study = names(expr)[k],
                surv.outcome = '${survival_type}',
                cancer.type = cancer_types[k],
                treatment = treatment_types[k]
            )
        } else if ('${survival_type}' == 'Response') {
            geneLogReg(
                dat.icb = expr[[k]],
                missing.perc = 0.5,
                const.int = 0.001,
                n.cutoff = 15,
                feature = "${params.gene_name}",
                study = names(expr)[k],
                n0.cutoff = 5, 
                n1.cutoff = 5,
                cancer.type = cancer_types[k],
                treatment = treatment_types[k]
            )
        }
    })

    # Combine results into one data frame
    assoc.res <- do.call(rbind, assoc.res)

    # Adjust p-values for multiple comparisons using the Benjamini-Hochberg method
    assoc.res\$FDR <- p.adjust(assoc.res\$Pval, method = "BH")

    # Meta-analysis for a gene across datasets
    if ('${survival_type}' == 'OS' || '${survival_type}' == 'PFS') {
        res_meta_pertreatment <- metaPerTreatmentfun(
            coef = assoc.res\$Coef, 
            se = assoc.res\$SE,
            study  = assoc.res\$Study, 
            pval = assoc.res\$Pval, 
            n = assoc.res\$N,
            cancer.type = assoc.res\$Cancer_type,
            treatment = assoc.res\$Treatment,
            feature = "${params.gene_name}", 
            treatment.spec = TRUE
        )
    } else if ('${survival_type}' == 'Response') {
        res_meta_pertreatment <- metaLogRegfun(
            coef = assoc.res\$Coef, 
            se = assoc.res\$SE,
            study  = assoc.res\$Study, 
            pval = assoc.res\$Pval, 
            n = assoc.res\$N,
            cancer.type = assoc.res\$Cancer_type,
            treatment = assoc.res\$Treatment,
            feature = "${params.gene_name}", 
            treatment.spec = TRUE
        )
    }

    res_meta_pertreatment <- data.frame(rbind(res_meta_pertreatment\$PD1\$meta_summery, res_meta_pertreatment\$CTLA4\$meta_summery))

    # Save the results to a CSV file
    write.csv(res_meta_pertreatment, file = "Meta_analysis_${survival_type}_${params.gene_name}_pertreatment.csv", row.names = FALSE)
    """
}


workflow {
    /*
    ===============================================================
    SECTION: Load RDA Data and Extract Expression and Clinical Data
    ===============================================================
    */

    // Set the specified study id and data directory
    icb_dat = file("${params.data_dir}/${params.study_id}.rda")

    // Extract expression and clinical data to CSV files
    extracted_data = LoadAndExtractData(icb_dat)

    // Define clinical and expression data files
    expr_file = extracted_data[0]
    clin_file = extracted_data[1]

    /*
    ========================================================
    SECTION: Gene Association Analysis
    ========================================================
    */
    
    /*
    --------------------------------------------------------
    A. SUBSECTION: OS or PFS
    --------------------------------------------------------
    */

    // Here survival_type is "OS". You can use "PFS" alternatively.
    // Also define your genes like genes = "CXCL9" or vector of genes as below

    genes = 'c("CXCL9", "CXCL10", "TIGIT", "CD83", "STAT1", "CXCL11", "CXCL13", "CD8A", "CTLA4")'
    survival_result_os = GeneAssociationSurvival(expr_file, clin_file, survival_type = "OS", genes)
    response_result = GeneAssociationResponse(expr_file, clin_file, genes)

    /*
    ========================================================
    SECTION: Meta Analysis - OS
    ========================================================
    */

    all_rdas = file(params.data_dir)

    // Perform meta-analysis for OS; you can use survival_type "PFS" or "Response"
    MetaAnalysis_pertreatment(all_rdas, survival_type = "OS")

}