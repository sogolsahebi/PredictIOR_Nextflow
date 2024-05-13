#!/usr/bin/env nextflow
nextflow.enable.dsl=2

params.study_id = 'ICB_Ravi__Lung__PD-L1'
params.data_dir = './data'
params.out_dir = './results'
params.rda_file = "${params.data_dir}/${params.study_id}.rda"

process LoadAndExtractData {
    container 'nextflow-r-env:latest'
    publishDir "./outputs", mode: 'copy'

    input:
    path rda_file

    output:
    tuple path("${params.out_dir}/${params.study_id}_expr.csv"), path("${params.out_dir}/${params.study_id}_clin.csv")

    script:
    """
    #! /usr/bin/env Rscript
    library(SummarizedExperiment)

    load("${rda_file}")

    expr <- assay(dat_icb)
    clin <- as.data.frame(colData(dat_icb))

    output_dir <- "${params.out_dir}"
    if (!dir.exists(output_dir)) {
        dir.create(output_dir, recursive = TRUE)
    }

    expr_file_path <- paste0(output_dir, "/", "${params.study_id}_expr.csv")
    clin_file_path <- paste0(output_dir, "/", "${params.study_id}_clin.csv")

    write.csv(expr, expr_file_path, row.names = TRUE)
    write.csv(clin, clin_file_path, row.names = FALSE)
    """
}

process ViewCSVContents {
    container 'nextflow-r-env:latest'

    input:
    tuple path(expr_file), path(clin_file)

    script:
    """
    #! /usr/bin/env Rscript
    library(data.table)

    expr_data <- fread("${expr_file}")
    clin_data <- fread("${clin_file}")

    """
}

workflow {
    icb_dat = file(params.rda_file)
    out = LoadAndExtractData(icb_dat)
    ViewCSVContents(out)
}
