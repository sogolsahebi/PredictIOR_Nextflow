#!/usr/bin/env nextflow
nextflow.enable.dsl=2

process RunRScript {
    output:
    path "*.rda"

    script:
    """
    Rscript ${baseDir}/R/getSummarizedExperiment.R
    """
}

workflow {
    RunRScript()
    output = Channel.fromPath("*.rda")
    output.collectFile(name: 'results.rda', storeDir: 'results')
}




// nextflow.enable.dsl=2

// process RunRScript {
//     output:
//     path "*.rda"

//     script:
//     """
//     Rscript /workspaces/predictIO-template/R/getSummarizedExperiment.R
//     """
// }

// workflow {
//     RunRScript()
//     output = Channel.fromPath("*.rda")
//     output.collectFile(name: 'results.rda', storeDir: '/workspaces/predictIO-template/results')
// }




// // Define script parameters
// params.dataDir = './data'
// params.outDir = './results'
// params.container = 'bioconductor/bioconductor_docker:RELEASE_3_12'

// log.info """
//     P R E D I C T I O - N F   P I P E L I N E
//     ===================================
//     dataDir       : ${params.dataDir}
//     outDir        : ${params.outDir}
//     container     : ${params.container}
// """

// // Define input data channel
// Channel
//     .fromPath("${params.dataDir}/*.rds")
//     .set { input_data_ch }

// process LoadData {
//     tag "Loading data"
//     publishDir "${params.outDir}", mode: 'copy', overwrite: true
//     container 'bioconductor/bioconductor_docker:RELEASE_3_12'
//     containerOptions = '-v /mnt/c/Users/sogol/predictIO-template/R:/R'

//     input:
//     path data_file

//     output:
//     path "${params.outDir}/${data_file.baseName}.rds"

//     script:
//     """
//     # Initialize Conda and install required packages
//     eval "$(conda shell.bash hook)"
//     conda install -c conda-forge -c bioconda -y bioconductor-multiassayexperiment

//     # Check for the presence of the R script
//     echo "Checking R script presence in mounted directory:"
//     ls -l /R/
//     Rscript --vanilla /R/getSummarizedExperiment.R $data_file ${data_file.baseName}.rds
//     """
// }



// // Main workflow that orchestrates the execution of the process
// workflow {
//     load_data = LoadData(input_data_ch)

//     load_data.view()
// }



// TODO: analaysis data.
// TODO: packages, redable, check : if the data csv file.
// Main workflow definition that orchestrates the execution of all processes




// // Process to perform meta analysis using getMetaAnalysis.R
// process PerformMetaAnalysis {
//     tag "Meta Analysis"
//     publishDir "${params.outDir}", mode: 'copy', overwrite: true
//     container params.container

//     input:
//     path data

//     output:
//     path 'meta_analysis_results.txt' into meta_analysis_results_ch

//     script:
//     """
//     Rscript --vanilla getMetaAnalysis.R $data meta_analysis_results.txt
//     """
// }

// // Process to calculate gene signature score using getGeneSigScore.R
// process CalculateGeneSigScore {
//     tag "Gene Signature Score"
//     publishDir "${params.outDir}", mode: 'copy', overwrite: true
//     container params.container

//     input:
//     path data

//     output:
//     path 'gene_sig_score_results.txt' into gene_sig_score_results_ch

//     script:
//     """
//     Rscript --vanilla getGeneSigScore.R $data gene_sig_score_results.txt
//     """
// }

// // Main workflow definition that orchestrates the execution of all processes
// workflow {
//     load_data = LoadData(input_data_ch)
//     meta_analysis = PerformMetaAnalysis(load_data.out)
//     gene_sig_score = CalculateGeneSigScore(load_data.out)
// }
