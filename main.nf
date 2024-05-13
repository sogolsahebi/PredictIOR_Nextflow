#!/usr/bin/env nextflow
nextflow.enable.dsl=2

params.study_id = 'ICB_Ravi__Lung__PD-L1'
params.data_dir = './data'
params.out_dir = './output'
params.rda_file = "${params.data_dir}/${params.study_id}.rda"

// Load RDA data, extract expression and clinical data, and write to CSV files.
process LoadAndExtractData {
    container 'nextflow-r-env:latest'
    publishDir "${params.out_dir}", mode: 'copy'

    input:
    path rda_file

    output:
    tuple path("${params.study_id}_expr.csv"), path("${params.study_id}_clin.csv")

    script:
    """
    #! /usr/bin/env Rscript
    library(SummarizedExperiment)

    load("${rda_file}")

    expr <- assay(dat_icb)
    clin <- as.data.frame(colData(dat_icb))

    write.csv(expr, "${params.study_id}_expr.csv", row.names = TRUE)
    write.csv(clin, "${params.study_id}_clin.csv", row.names = FALSE)
    """
}




workflow {
    // load rda file and extact clin and expr
    icb_dat = file(params.rda_file)
    LoadAndExtractData(icb_dat)

    //


}
