#!/usr/bin/env nextflow
nextflow.enable.dsl=2

params.study_id = 'ICB_Ravi__Lung__PD-L1'
params.data_dir = './data'
params.out_dir = './output'
params.rda_file = "${params.data_dir}/${params.study_id}.rda"

//single gene - vector , feautre input vs all genes nrow(expr)

// Process to load RDA data and extract expression and clinical data
process LoadAndExtractData {
    container 'nextflow-rmd-env'
    publishDir "${params.out_dir}", mode: 'copy'

    input:
    path rda_file

    output:
    path "${params.study_id}_expr.csv"
    path "${params.study_id}_clin.csv"

    script:
    """
    #!/usr/bin/env Rscript
    library(SummarizedExperiment)

    load("${rda_file}")

    expr <- assay(dat_icb)
    clin <- as.data.frame(colData(dat_icb))

    write.csv(expr, "${params.study_id}_expr.csv", row.names = TRUE)
    write.csv(clin, "${params.study_id}_clin.csv", row.names = FALSE)
    """
}

// Process to run a specific R script for gene association analysis
//1.Gene Association  with 'Os'
process GeneAssociationOs {
    tag "${params.study_id}"
    container 'nextflow-rmd-env'
    publishDir "${params.out_dir}", mode: 'copy'

    input:
    path expr_file
    path clin_file

    output:
    path "${params.study_id}_cox_os_results.csv"

    script:
    def studyParts = params.study_id.split("__")
    """
    #!/usr/bin/env Rscript
    library(SummarizedExperiment)

    args <- commandArgs(trailingOnly = TRUE)
    study_id_parts <- unlist(strsplit(args[1], " "))

    expr <- read.csv("${expr_file}", row.names = 1)
    clin <- read.csv("${clin_file}")

    source('/R_scripts/getGeneAssociation.R')
    source('/R_scripts/getHR.R')

    # TODO: gene_level , adding Signiture level?
    # for single gene?

    results <- lapply(1:100, function(i) {
      geneSurvCont(
        dat.icb = expr,
        clin = clin,
        time.censor = 36,
        missing.perc = 0.5,
        const.int = 0.001,
        n.cutoff = 15,
        feature = rownames(expr)[i],  # TODO: single gene - all genes.          
        study = paste(study_id_parts[0], study_id_parts[1], study_id_parts[2], sep='__'),
        surv.outcome = 'OS',
        cancer.type = study_id_parts[1],
        treatment = study_id_parts[2]
      )
    })

    # Convert list to a data frame
    results_df <- do.call(rbind, results)

    write.csv(results_df, file = "${params.study_id}_cox_os_results.csv", row.names = FALSE)
    """
}

//2.Gene Association  with Immunotherapy response (R vs NR).
process LogisticRegression {
    tag "${params.study_id}"
    container 'nextflow-rmd-env'
    publishDir "${params.out_dir}", mode: 'copy'

    input:
    path expr_file
    path clin_file

    output:
    path "${params.study_id}_results_logreg.csv"

    script:
    def studyParts = params.study_id.split("__")
    """
    #!/usr/bin/env Rscript
    library(SummarizedExperiment)

    args <- commandArgs(trailingOnly = TRUE)
    study_id_parts <- unlist(strsplit(args[1], "__"))

    expr <- read.csv("${expr_file}", row.names = 1)
    clin <- read.csv("${clin_file}")

    source('/R_scripts/getGeneAssociation.R')
    source('/R_scripts/getHR.R')

    logreg <- lapply(1:min(200, nrow(expr)), function(i) {
      res <- geneLogReg(
        dat.icb = expr,
        clin = clin,
        missing.perc = 0.5,
        const.int = 0.001,
        n.cutoff = 15,
        feature = rownames(expr)[i],
        study = paste(study_id_parts[0], study_id_parts[1], study_id_parts[2], sep = '__'),
        n0.cutoff = 3,
        n1.cutoff = 3,
        cancer.type = study_id_parts[1],
        treatment = study_id_parts[2]
      )
      res
    })

    # Convert list to a data frame
    logreg_df <- do.call(rbind, logreg)

    # Save the results to a CSV file
    write.csv(logreg_df, file = "${params.study_id}_results_logreg.csv", row.names = FALSE)
    """
}

//3. PFS

workflow {
    // Entry point for the RDA data file
    icb_dat = file(params.rda_file)

    // Load RDA data and extract expression and clinical data to CSV files
    extracted_data = LoadAndExtractData(icb_dat)

    // Use the extracted CSV files for gene association analysis
    GeneAssociationOs(
        extracted_data[0],  // expr.csv
        extracted_data[1]   // clin.csv
    )

    // Use the extracted CSV files for logistic regression analysis
    LogisticRegression(extracted_data[0], extracted_data[1])
}