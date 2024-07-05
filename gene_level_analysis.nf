#!/usr/bin/env nextflow
nextflow.enable.dsl=2

params.study_id = 'ICB_small_Mariathasan'
params.data_dir = './ICB_data'
params.out_dir = './output'
params.rda_file = "${params.data_dir}/${params.study_id}.rda"
cancer_type = 'Bladder'
treatment = 'PD-1/PD-L1'

/*
========================================================
SECTION: Load Immunotherapy Datasets
========================================================
*/

/*
Load public clinical multimodal immunotherapy datasets from GitHub and ORCESTRA for transparent biomarker discovery in immunotherapy response. For RNA profiles, we use log2-transformed TPM data from protein-coding genes, filtering out genes with zero expression in at least 50% of samples. Only studies with at least 20 patients are included.

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
    library(SummarizedExperiment)

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
======================================================================================
SECTION: Biomarkers and Immunotherapy Response Association - Gene Association Analysis
======================================================================================
*/

/*
Assessing the association of specific biomarkers with:
- Immunotherapy response (R vs NR)
- Immunotherapy survival (for time-to-event analysis: OS and PFS)

P-values are corrected for multiple testing using the Benjamini-Hochberg (FDR) method, with significance set at p-values or FDR ≤ 5%.

Note: Coef represents the log hazard ratio (logHR) or log odds ratio (logOR) depending on the analysis.
*/

/*
--------------------------------------------------------
SUBSECTION: Overall survival (OS)
--------------------------------------------------------
*/

genes = 'c("CXCL9", "CXCL10", "TIGIT", "CD83", "STAT1", "CXCL11", "CXCL13", "CD8A", "CTLA4")'

process GeneAssociationOS {
    tag "${params.study_id}"
    container 'nextflow-env:latest'
    publishDir "${params.out_dir}", mode: 'copy'

    input:
    path expr_file
    path clin_file

    output:
    path "${params.study_id}_cox_os_genes.csv"

    script:
    """
    #!/usr/bin/env Rscript
    library(SummarizedExperiment)
    library(PredictioR)
    library(survival)
    library(BiocGenerics)
    library(S4Vectors)
    library(IRanges)
    library(GenomeInfoDb)
    library(Biobase)
    
    expr <- read.csv("${expr_file}", row.names = 1)
    clin <- as.data.frame(colData(${params.study_id}))

    cox_os <- geneSurvCont(
        dat.icb = expr,
        clin = clin,
        time.censor = 36,
        missing.perc = 0.5,
        const.int = 0.001,
        n.cutoff = 15,
        feature = ${genes},
        study = "${params.study_id}",
        surv.outcome = 'OS',
        cancer.type = "${cancer_type}",
        treatment = "${treatment}"
    )

    # Adjust p-values for multiple testing
    cox_os\$FDR <- p.adjust(cox_os\$Pval, method = "BH")
    cox_os <- cox_os[order(cox_os\$FDR), ]
    write.csv(cox_os, file = "${params.study_id}_cox_os_genes.csv", row.names = FALSE)
    """
}

/*
--------------------------------------------------------
SUBSECTION: Progression-free survival (PFS)
--------------------------------------------------------
*/

process GeneAssociationPFS {
    tag "${params.study_id}"
    container 'nextflow-env:latest'
    publishDir "${params.out_dir}", mode: 'copy'

    input:
    path expr_file
    path clin_file

    output:
    path "${params.study_id}_cox_pfs_genes.csv"

    script:
    """
    #!/usr/bin/env Rscript
    library(SummarizedExperiment)
    library(PredictioR)
    library(survival)
    library(BiocGenerics)
    library(S4Vectors)
    library(IRanges)
    library(GenomeInfoDb)
    library(Biobase)
    
    expr <- read.csv("${expr_file}", row.names = 1)
    clin <- as.data.frame(colData(${params.study_id}))

    cox_pfs <- geneSurvCont(
        dat.icb = expr,
        clin = clin,
        time.censor = 36,
        missing.perc = 0.5,
        const.int = 0.001,
        n.cutoff = 15,
        feature = ${genes},
        study = "${params.study_id}",
        surv.outcome = 'PFS',
        cancer.type = "${cancer_type}",
        treatment = "${treatment}"
    )

    # Adjust p-values for multiple testing
    cox_pfs\$FDR <- p.adjust(cox_pfs\$Pval, method = "BH")
    cox_pfs <- cox_pfs[order(cox_pfs\$FDR), ]
    write.csv(cox_pfs, file = "${params.study_id}_cox_pfs_genes.csv", row.names = FALSE)
    """
}

/*
--------------------------------------------------------
SUBSECTION: Immunotherapy response (R vs NR)
--------------------------------------------------------
*/

process GeneAssociationResponse_RvsNR {
    tag "${params.study_id}"
    container 'nextflow-env:latest'
    publishDir "${params.out_dir}", mode: 'copy'

    input:
    path expr_file
    path clin_file

    output:
    path "${params.study_id}_Response_RvsNR.csv"

    script:
    """
    #!/usr/bin/env Rscript
    library(SummarizedExperiment)
    library(PredictioR)

    expr <- read.csv("${expr_file}", row.names = 1)
    clin <- read.csv("${clin_file}")
    
    logreg <- geneLogReg(
        dat.icb = expr,
        clin = clin,
        missing.perc = 0.5,
        const.int = 0.001,
        n.cutoff = 15,
        feature = c("CXCL9", "CXCL10", "TIGIT", "CD83", "STAT1", "CXCL11", "CXCL13", "CD8A", "CTLA4"),
        study = "${params.study_id}", 
        n0.cutoff = 5, 
        n1.cutoff = 5,
        cancer.type = "${cancer_type}",
        treatment = "${treatment}"
    )

    # Adjust P-values and sort by FDR
    logreg <- logreg[order(logreg\$FDR <- p.adjust(logreg\$Pval, method = "BH"))]

    # Save as CSV file
    write.csv(logreg, file = "${params.study_id}_Response_RvsNR.csv", row.names = FALSE)
    """
}

/*
========================================================
SECTION: Meta Analysis Section
========================================================
*/

/*
The following clinical multimodal immunotherapy datasets are publicly available on the Github. These datasets are used in biomarker discovery for immunotherapy response through pan-cancer, cancer-specific, and treatment-specific analyses

Links:
- GitHub: https://github.com/bhklab/PredictioR/tree/main/data
*/

// Set the cancer type and treatment for each dataset respectively.
params.cancer_types = '["Melanoma", "Pancreas", "Melanoma", "Bladder", "Melanoma", "Melanoma", "Kidney", "Melanoma"]'
params.treatment_types = '["PD-1/PD-L1", "PD-1/PD-L1", "PD-1/PD-L1", "PD-1/PD-L1", "CTLA4", "IO+combo", "PD-1/PD-L1", "CTLA4"]'
params.gene_name = "CXCL9"

/*
-----------------------------------------------------------------------
SUBSECTION: Aggregating Associations through Meta-analysis (Pan-cancer)
-----------------------------------------------------------------------
Set the directory that contains all the .rda files for meta-analysis. Previously defined as params.data_dir = './data'. There are 8 .rda files: "ICB_Liu", "ICB_Padron", "ICB_Hugo", "ICB_Mariathasan", "ICB_Nathanson", "ICB_Riaz", "ICB_Miao", "ICB

_Van_Allen"

As an example, for the gene CXCL9, to generalize the association with immunotherapy survival, we apply a meta-analysis approach to integrate findings across datasets for pan-cancer analysis.
*/

process MetaAnalysisOS_PanCancer {
    tag {params.gene_name }
    container 'nextflow-env:latest'
    publishDir "${params.out_dir}", mode: 'copy'

    input:
    path rda_files_path

    output:
    path "Meta_analysis_os_${params.gene_name}_association.csv"
    path "Meta_analysis_os_${params.gene_name}_pancancer.csv"

    script:
    """
    #!/usr/bin/env Rscript
    library(SummarizedExperiment)
    library(PredictioR)
    library(survival)
    library(meta)
    library(jsonlite) 

    # Load all .rda files 
    lapply(list.files(path = '${rda_files_path}', pattern = '*.rda', full.names = TRUE), function(file) {
        load(file, envir = .GlobalEnv)
    })

    # Create a list of the loaded objects using their actual names
    loaded_objects <- ls(pattern = "^ICB_small_")
    expr <- mget(loaded_objects, envir = .GlobalEnv)

    # Define the cancer types and treatment types vectors
    cancer_types <- fromJSON('${params.cancer_types}')
    treatment_types <- fromJSON('${params.treatment_types}')

    # Apply a function over the loaded datasets to perform survival analysis
    assoc.res <- lapply(1:length(expr), function(k){
        geneSurvCont(
            dat.icb = expr[[k]],
            time.censor = 36,
            missing.perc = 0.5,
            const.int = 0.001,
            n.cutoff = 15,
            feature = "${params.gene_name}",
            study = names(expr)[k],
            surv.outcome = 'OS',
            cancer.type = cancer_types[k],
            treatment = treatment_types[k]
        )
    })

    # Combine results into one data frame
    assoc.res <- do.call(rbind, assoc.res)

    # Adjust p-values for multiple comparisons using the Benjamini-Hochberg method
    assoc.res\$FDR <- p.adjust(assoc.res\$Pval, method = "BH")

    # Meta-analysis for a gene across datasets
    res_meta_pancancer <- metafun(
        coef = assoc.res\$Coef, 
        se = assoc.res\$SE,
        study  = assoc.res\$Study, 
        pval = assoc.res\$Pval, 
        n = assoc.res\$N,
        cancer.type = assoc.res\$Cancer_type,
        treatment = assoc.res\$Treatment,
        feature = "${params.gene_name}", 
        cancer.spec = FALSE, 
        treatment.spec = FALSE
    )
        
    # Meta-analysis results
    res_meta_pancancer <- data.frame(res_meta_pancancer)

    # Save the results to a CSV file
    write.csv(assoc.res, file = "Meta_analysis_os_${params.gene_name}_association.csv", row.names = FALSE)
    write.csv(res_meta_pancancer, file = "Meta_analysis_os_${params.gene_name}_pancancer.csv", row.names = FALSE)
    """
}

/*
-----------------------------------------------------------------------
SUBSECTION: Aggregating Associations through Meta-analysis (Per-cancer)
-----------------------------------------------------------------------

For cancer-specific analysis, consider meta-analysis when there are at least 3 datasets.
*/

process MetaAnalysisOS_PerCancer {
    tag {params.gene_name }
    container 'nextflow-env:latest'
    publishDir "${params.out_dir}", mode: 'copy'

    input:
    path rda_files_path

    output:
    path "Meta_analysis_os_${params.gene_name}_percancer.csv"

    script:
    """
    #!/usr/bin/env Rscript
    library(SummarizedExperiment)
    library(PredictioR)
    library(survival)
    library(meta)
    library(jsonlite)

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

    # Apply a function over the loaded datasets to perform survival analysis
    assoc.res <- lapply(1:length(expr), function(k){
        geneSurvCont(
            dat.icb = expr[[k]],
            time.censor = 36,
            missing.perc = 0.5,
            const.int = 0.001,
            n.cutoff = 15,
            feature = "${params.gene_name}",
            study = names(expr)[k],
            surv.outcome = 'OS',
            cancer.type = cancer_types[k],
            treatment = treatment_types[k]
        )
    })

    # Combine results into one data frame
    assoc.res <- do.call(rbind, assoc.res)

    # Adjust p-values for multiple comparisons using the Benjamini-Hochberg method
    assoc.res\$FDR <- p.adjust(assoc.res\$Pval, method = "BH")

    # Treatment-specific meta-analysis for a gene across datasets
    res_meta_percancer <- metaPerCanfun(
        coef = assoc.res\$Coef, 
        se = assoc.res\$SE,
        study  = assoc.res\$Study, 
        pval = assoc.res\$Pval, 
        n = assoc.res\$N,
        cancer.type = assoc.res\$Cancer_type,
        treatment = assoc.res\$Treatment,
        feature = "${params.gene_name}", 
        cancer.spec = TRUE
    )

    res_meta_percancer <- data.frame(rbind(res_meta_percancer\$Melanoma\$meta_summery, res_meta_percancer\$Other\$meta_summery))

    # Save the results to a CSV file
    write.csv(res_meta_percancer, file = "Meta_analysis_os_${params.gene_name}_percancer.csv", row.names = FALSE)
    """
}

/*
--------------------------------------------------------
SUBSECTION: Aggregating Associations through Meta-analysis (Per-treatment)
--------------------------------------------------------

As an example, for the gene CXCL9, to generalize the association with immunotherapy survival, we apply a meta-analysis approach to integrate findings across datasets for pan-cancer analysis.
*/

process MetaAnalysisOS_PerTreatment {
    tag {params.gene_name}
    container 'nextflow-env:latest'
    publishDir "${params.out_dir}", mode: 'copy'

    input:
    path rda_files_path

    output:
    path "Meta_analysis_os_${params.gene_name}_pertreatment.csv"

    script:
    """
    #!/usr/bin/env Rscript
    library(SummarizedExperiment)
    library(PredictioR)
    library(survival)
    library(meta)
    library(jsonlite)

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

    # Apply a function over the loaded datasets to perform survival analysis
    assoc.res <- lapply(1:length(expr), function(k){
        geneSurvCont(
            dat.icb = expr[[k]],
            time.censor = 36,
            missing.perc = 0.5,
            const.int = 0.001,
            n.cutoff = 15,
            feature = "${params.gene_name}",
            study = names(expr)[k],
            surv.outcome = 'OS',
            cancer.type = cancer_types[k],
            treatment = treatment_types[k]
        )
    })

    # Combine results into one data frame
    assoc.res <- do.call(rbind, assoc.res)

    # Adjust p-values for multiple comparisons using the Benjamini-Hochberg method
    assoc.res\$FDR <- p.adjust(assoc.res\$Pval, method = "BH")

    # Treatment-specific meta-analysis for a gene across datasets
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
        
    # Summary of meta-analysis results for each treatment type
    res_meta_pertreatment <- data.frame(rbind(res_meta_pertreatment\$`PD-1/PD-L1`\$meta_summery, res_meta_pertreatment\$Other\$meta_summery))

    # Save the results to a CSV file
    write.csv(res_meta_pertreatment, file = "Meta_analysis_os_${params.gene_name}_pertreatment.csv", row.names = FALSE)
    """
}


workflow {
    /*
    ===============================================================
    SECTION: Load RDA Data and Extract Expression and Clinical Data
    ===============================================================
    */

    // Set the specified study id and data directory
    icb_dat = file(params.rda_file)

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

    // Perform gene association analysis

    /*
    --------------------------------------------------------
    A. SUBSECTION: OS and PFS
    --------------------------------------------------------
    */

    GeneAssociationOS(expr_file, clin_file)
    GeneAssociationPFS(expr_file, clin_file)

    /*
    --------------------------------------------------------
    B. SUBSECTION: Immunotherapy response (R vs NR)
    --------------------------------------------------------
    */

    GeneAssociationResponse_RvsNR(expr_file, clin_file)

    /*
    ========================================================
    SECTION: Meta Analysis - OS
    ========================================================
    */

    all_rdas = file(params.data_dir)

    // Perform Meta Analysis (OS) using loaded RDA data
    MetaAnalysisOS_PanCancer(all_rdas)
    MetaAnalysisOS_PerCancer(all_rdas)
    MetaAnalysisOS_PerTreatment(all_rdas)
}
