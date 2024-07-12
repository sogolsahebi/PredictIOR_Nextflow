#!/usr/bin/env nextflow
nextflow.enable.dsl=2

params.study_id = 'ICB_small_Padron' // Set study_id (name of the .rda file in ./ICB_data directory)
params.data_dir = './ICB_data'
params.out_dir = './output/gene_level_output'
params.cancer_type = 'Pancreas' // Set the cancer type of the study you chose
params.treatment = 'PD-1/PD-L1' // Set the treatment of the study you chose


/*
Note
Another input example would be:

params.study_id = 'ICB_small_Liu' 
params.cancer_type = 'Melanoma' 
params.treatment = 'PD-1/PD-L1'
/*

log.info """
P R E D I C T I O - N F   P I P E L I N E (Gene Level Analysis)
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
Public clinical multimodal immunotherapy datasets are from GitHub or ORCESTRA for 
transparent biomarker discovery in immunotherapy response. For RNA profiles, we use
log2-transformed TPM data from protein-coding genes, filtering out genes with zero 
expression in at least 50% of samples. Only studies with at least 20 patients are included.

For example, the Padron dataset (PMID: 35662283) includes RNA expression, clinical data,
and gene metadata for 45 patients with 18,459 protein-coding genes, focused on Pancreas 
cancer and PD-1/PD-L1 treatment. Here, a subset of the dataset used and shared (ICB_small_Padron).

Links:
- GitHub: https://github.com/bhklab/PredictioR/tree/main/data
- ORCESTRA: https://www.orcestra.ca/clinical_icb

The detailed clinical characteristics of the Mariathasan dataset include cancer type, age, sex, response (R vs. NR), overall survival (OS) and Progression-free survival (PFS).
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
Notes:
For gene-level analysis in the R script, you can use dat.icb in two ways:

Option 1: Use dat.icb = expr (data frame) with clin = clin (data frame) for clinical data.
Option 2: Load the RDA file with load("${rda_file}"), then use dat.icb = '${params.study_id}' (SummarizedExperiment object) with clin = NULL.
*/

/*
========================================================
SECTION: Biomarkers and Immunotherapy Response Association 
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
----------------
SUBSECTION: OS
----------------
*/

process GeneAssociationOS {
    tag "${params.study_id}"
    container 'nextflow-env:latest'
    publishDir "${params.out_dir}", mode: 'copy'

    input:
    path expr_file
    path clin_file

    val genes

    output:
    path("${params.study_id}_cox_os.csv")

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
        surv.outcome = "OS",
        cancer.type = "${params.cancer_type}",
        treatment = "${params.treatment}"
    )

    # Adjust p-values for multiple testing
    cox_result\$FDR <- p.adjust(cox_result\$Pval, method = "BH")
    cox_result <- cox_result[order(cox_result\$FDR), ]
    write.csv(cox_result, file = "${params.study_id}_cox_os.csv", row.names = FALSE)
    """
}

/*
-----------------------
SUBSECTION: PFS
----------------------
*/

process GeneAssociationPFS {
    tag "${params.study_id}"
    container 'nextflow-env:latest'
    publishDir "${params.out_dir}", mode: 'copy'

    input:
    path expr_file
    path clin_file
    val genes

    output:
    path("${params.study_id}_cox_pfs.csv")

    script:
    """
    #!/usr/bin/env Rscript
    source('/R/load_libraries.R')
    
    expr <- read.csv("${expr_file}", row.names = 1)
    clin <- read.csv("${clin_file}")

    cox_result <- geneSurvCont(
        dat.icb = expr,
        clin = clin,
        time.censor = 24,
        missing.perc = 0.5,
        const.int = 0.001,
        n.cutoff = 15,
        feature = ${genes},
        study = "${params.study_id}",
        surv.outcome = "PFS",
        cancer.type = "${params.cancer_type}",
        treatment = "${params.treatment}"
    )

    # Adjust p-values for multiple testing
    cox_result\$FDR <- p.adjust(cox_result\$Pval, method = "BH")
    cox_result <- cox_result[order(cox_result\$FDR), ]
    write.csv(cox_result, file = "${params.study_id}_cox_pfs.csv", row.names = FALSE)
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
        n0.cutoff = 3, 
        n1.cutoff = 3,
        cancer.type = "${params.cancer_type}",
        treatment = "${params.treatment}"
    )

    # Adjust P-values and sort by FDR
    logreg <- logreg[order(logreg\$FDR <- p.adjust(logreg\$Pval, method = "BH"))]

    # Save as CSV file
    write.csv(logreg, file = "${params.study_id}_Response.csv", row.names = FALSE)
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
    A. SUBSECTION: OS 
    --------------------------------------------------------
    */

    // You can choose your genes as a single gene or a vector of genes as shown below
    genes = 'c("CXCL9", "CXCL10", "TIGIT", "CD83", "STAT1", "CXCL11", "CXCL13", "CD8A", "CTLA4")'
    survival_result_os = GeneAssociationOS(expr_file, clin_file, genes)

    // If using a single gene, uncomment the following lines
    //gene = "CXCL9"
    //survival_result_os = GeneAssociationOS(expr_file, clin_file, gene)

    /*
    --------------------------------------------------------
    B. SUBSECTION: PFS 
    --------------------------------------------------------
    */

    survival_result_os = GeneAssociationPFS(expr_file, clin_file, genes)

    /*
    --------------------------------------------------------
    C. SUBSECTION: Immunotherapy response (R vs NR) 
    --------------------------------------------------------
    */
    response_result = GeneAssociationResponse(expr_file, clin_file, genes)


}
