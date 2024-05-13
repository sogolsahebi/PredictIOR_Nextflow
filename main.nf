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

    output:
    file 'gene_association_output.txt'

    script:
    """
    #!/usr/bin/env Rscript
    if (file.exists("R/getGeneAssociation.R")) {
        source("R/getGeneAssociation.R")
        message("getGeneAssociation.R sourced successfully.")
    } else {
        message("Error: R/getGeneAssociation.R not found.")
    }
    """
}



workflow {

    // Define study parts for use in analysis within the workflow block
    study_id_parts = params.study_id.split("__")

    // Entry point for the RDA data file
    icb_dat = file(params.rda_file)

    // Load RDA data and extract expression and clinical data to CSV files
    LoadAndExtractData(icb_dat)

    // Use extracted data to analyze gene survival associations
    GeneAssociationOs(LoadAndExtractData.out, params.study_id, study_id_parts)
}

