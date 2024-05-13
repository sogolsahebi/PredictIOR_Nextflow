#!/usr/bin/env nextflow
nextflow.enable.dsl=2

params.study_id = 'ICB_Ravi__Lung__PD-L1'
params.data_dir = './data'
params.out_dir = './output'
params.rda_file = "${params.data_dir}/${params.study_id}.rda"

process LoadAndExtractData {
    container 'nextflow-r-env:latest'
    publishDir "${params.out_dir}", mode: 'copy'

    input:
    path rda_file

    output:
    tuple path("${params.study_id}_expr.csv"), path("${params.study_id}_clin.csv")

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

process GeneAssociationOs {
    container 'nextflow-r-env:latest'
    publishDir "${params.out_dir}", mode: 'copy'

    input:
    tuple path(expr_csv), path(clin_csv)
    val study_id

    output:
    path("cox_os_results.csv")

    script:
    """
    #!/usr/bin/env Rscript
    expr <- read.csv("${expr_csv}", row.names = 1)
    clin <- read.csv("${clin_csv}")

    # Check if the source file can be accessed and the function can be sourced without error
    box::use("R/getGeneAssociation.R[...]")


    # Ensure the split and paste are working as expected
    study_id_parts <- unlist(strsplit("${study_id}", "__"))
    print(study_id_parts)

    # Placeholder for actual function call
    # Insert actual function call here once above checks are confirmed to be correct
    """
}

workflow {

    // Define study parts for use in analysis
    study_id_parts = params.study_id.split("__")

    // Entry point for the RDA data file
    icb_dat = file(params.rda_file)

    // Load RDA data and extract expression and clinical data to CSV files
    LoadAndExtractData(icb_dat)

    // Use extracted data to analyze gene survival associations
    GeneAssociationOs(LoadAndExtractData.out, params.study_id)

}
