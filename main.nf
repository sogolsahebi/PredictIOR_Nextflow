#!/usr/bin/env nextflow
nextflow.enable.dsl=2

params.study_id = 'ICB_Ravi__Lung__PD-L1'
params.data_dir = './data'
params.out_dir = './output'
params.rda_file = "${params.data_dir}/${params.study_id}.rda"

// Process to load RDA data and extract expression and clinical data
process LoadAndExtractData {
    container 'sogolsahebi/nextflow-r-env:latest'
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
process GeneAssociationOs {
    tag "${rda_file.baseName}"
    container 'sogolsahebi/nextflow-r-env:latest'
    publishDir "${params.out_dir}", mode: 'copy'

    input:
    path rda_file

    output:
    path "${rda_file.baseName}_results.csv"

    script:
    def studyParts = params.study_id.split("__")
    """
    #!/usr/bin/env Rscript
    library(SummarizedExperiment)

    args <- commandArgs(trailingOnly = TRUE)
    study_id_parts <- unlist(strsplit(args[1], " "))

    load("${rda_file}")
    expr <- assay(dat_icb)
    clin <- as.data.frame(colData(dat_icb))

    source('/R_scripts/getGeneAssociation.R')
    source('/R_scripts/getHR.R')

    results <- lapply(1:100, function(i) {
      geneSurvCont(
        dat.icb = expr,
        clin = clin,
        time.censor = 36,
        missing.perc = 0.5,
        const.int = 0.001,
        n.cutoff = 15,
        feature = rownames(expr)[i],
        study = paste(study_id_parts[1], study_id_parts[2], study_id_parts[3], sep='__'),
        surv.outcome = 'OS',
        cancer.type = study_id_parts[2],
        treatment = study_id_parts[3]
      )
    })

    # Convert list to a data frame
    results_df <- do.call(rbind, results)  # Make sure to use 'results' not 'results_list'
    write.csv(results_df, file = 'ICB_Ravi__Lung__PD-L1_results.csv', row.names = FALSE)
    """
}


workflow {
    study_id_parts = params.study_id.split("__")

    // Entry point for the RDA data file
    icb_dat = file(params.rda_file)

    // Load RDA data and extract expression and clinical data to CSV files
    LoadAndExtractData(icb_dat)
    GeneAssociationOs(icb_dat)
}

