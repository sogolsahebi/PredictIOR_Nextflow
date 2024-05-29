#!/usr/bin/env nextflow
nextflow.enable.dsl=2

params.study_id = 'ICB_Ravi__Lung__PD-L1'
params.data_dir = './data'
params.out_dir = './output'
params.rda_file = "${params.data_dir}/${params.study_id}.rda"

// Extract cancer type and treatment type of this dataset
cancer_type = params.study_id.split('__')[1]
treatment_type = params.study_id.split('__')[2]
println("Extracted Cancer Type: ${cancer_type}, Treatment Type: ${treatment_type} for ${params.study_id} ")

// Process to load RDA data and extract expression and clinical data
process LoadAndExtractData {
    container 'sogolsahebi/nextflow-rmd-env'
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

////////// Gene Association Anlaysis //////////

params.gene_nums = 100  // Default number of genes to process

//1.Gene Association with 'OS'
process GeneAssociationOS {
    tag "${params.study_id}"
    container 'sogolsahebi/nextflow-rmd-env'
    publishDir "${params.out_dir}", mode: 'copy'

    input:
    path expr_file
    path clin_file

    output:
    path "${params.study_id}_cox_os_results.csv"

    script:
    """
    #!/usr/bin/env Rscript
    library(SummarizedExperiment)

    expr <- read.csv("${expr_file}", row.names = 1)
    clin <- read.csv("${clin_file}")

    source('/R_scripts/getGeneAssociation.R')
    source('/R_scripts/getHR.R')

    cox_os <- lapply(1:${params.gene_nums}, function(i) {
        geneSurvCont(
            dat.icb = expr,
            clin = clin,
            time.censor = 36,
            missing.perc = 0.5,
            const.int = 0.001,
            n.cutoff = 15,
            feature = rownames(expr)[i],
            study = '${params.study_id}',
            surv.outcome = 'OS',
            cancer.type = '${cancer_type}',
            treatment = '${treatment_type}'
        )
    })

    cox_os <- do.call(rbind, cox_os)
    # Additional filtering and adjustments as needed
    cox_os <- cox_os[!is.na(cox_os\$Gene), ]
    cox_os\$FDR <- p.adjust(cox_os\$Pval, method = "BH")

    write.csv(cox_os, file = "${params.study_id}_cox_os_results.csv", row.names = FALSE)
    """
}

//2.Gene Association with 'PFS'
process GeneAssociationPFS {
    tag "${params.study_id}"
    container 'sogolsahebi/nextflow-rmd-env'
    publishDir "${params.out_dir}", mode: 'copy'

    input:
    path expr_file
    path clin_file

    output:
    path "${params.study_id}_cox_pfs_results.csv"

    script:
    """
    #!/usr/bin/env Rscript
    library(SummarizedExperiment)

    expr <- read.csv("${expr_file}", row.names = 1)
    clin <- read.csv("${clin_file}")

    source('/R_scripts/getGeneAssociation.R')
    source('/R_scripts/getHR.R')

    cox_pfs <- lapply(1:${params.gene_nums}, function(i) {
      geneSurvCont(
        dat.icb = expr,
        clin = clin,
        time.censor = 36,
        missing.perc = 0.5,
        const.int = 0.001,
        n.cutoff = 15,
        feature = rownames(expr)[i],
        study = '${params.study_id}',
        surv.outcome = 'PFS',
        cancer.type = '${cancer_type}',
        treatment = '${treatment_type}'
      )
    })

    cox_pfs <- do.call(rbind, cox_pfs)
    
    # Additional filtering and adjustments as needed
    cox_pfs <- cox_pfs[!is.na(cox_pfs\$Gene), ]
    cox_pfs\$FDR <- p.adjust(cox_pfs\$Pval, method = "BH")

    write.csv(cox_pfs, file = "${params.study_id}_cox_pfs_results.csv", row.names = FALSE)
    """
}

//2.Gene Association  with Immunotherapy response (R vs NR).
process LogisticRegression {
    tag "${params.study_id}"
    container 'sogolsahebi/nextflow-rmd-env'
    publishDir "${params.out_dir}", mode: 'copy'

    input:
    path expr_file
    path clin_file

    output:
    path "${params.study_id}_logreg_results.csv"

    script:
    """
    #!/usr/bin/env Rscript
    library(SummarizedExperiment)

    expr <- read.csv("${expr_file}", row.names = 1)
    clin <- read.csv("${clin_file}")

    source('/R_scripts/getGeneAssociation.R')
    source('/R_scripts/getHR.R')

    logreg <- lapply(1:${params.gene_nums}, function(i) {
      res <- geneLogReg(
        dat.icb = expr,
        clin = clin,
        missing.perc = 0.5,
        const.int = 0.001,
        n.cutoff = 15,
        feature = rownames(expr)[i],
        study = '${params.study_id}',
        n0.cutoff = 3,
        n1.cutoff = 3,
        cancer.type = '${cancer_type}',
        treatment = '${treatment_type}'
      )

      # Filter out the results that do not have a coefficient (Coef) value.
      res <- res[!is.na(res\$Coef), ]
      res
    })

    # Convert list to a data frame
    logreg_df <- do.call(rbind, logreg)

    # Save the results to a CSV file
    write.csv(logreg_df, file = "${params.study_id}_logreg_results.csv", row.names = FALSE)
    """
}

// Aggregating Associations (OS) through Meta-analysis (Pan-cancer)
params.rda_files_dir = "/workspace/PredictIOR_Nextflow/data"
params.gene_name = "CXCL9"


////////// Meta Analysis Section //////////

// Process to load and prepare RDA files
process LoadAllData{
    container 'sogolsahebi/nextflow-rmd-env'
    publishDir "${params.out_dir}", mode: 'copy'

    input:
    path rda_files_path

    output:
    path "prepared_data.rds"

    script:
    """
    #!/usr/bin/env Rscript

    # Load each RDA file and extract data
    list_rda <- lapply(list.files(path = '${rda_files_path}', pattern = '*.rda', full.names = TRUE), function(file) {
        load(file)
        dat_icb
    })

    # Assign study names to the list elements
    study_names <- substr(list.files('${rda_files_path}'), 5, nchar(list.files('${rda_files_path}')) - 4)
    names(list_rda) <- study_names

    # Save the list_rda and study_names for further processing
    saveRDS(list(list_rda, study_names), file = "prepared_data.rds")
    """
}

process MetaAnalysisOS{
    container 'sogolsahebi/nextflow-rmd-env'
    publishDir "${params.out_dir}", mode: 'copy'

    input:
    path prepared_data

    output:
    path "${params.study_id}_meta_analysis_os.csv"

    script:
    """
    #!/usr/bin/env Rscript

    # Load prepared data
    prepared_data <- readRDS('${prepared_data}')
    list_rda <- prepared_data[[1]]
    study_names <- prepared_data[[2]]

    # Calling the functions
    source('/R_scripts/getGeneAssociation.R')
    source('/R_scripts/getHR.R')

    # Apply a function over the loaded datasets to perform survival analysis
    cox_os <- lapply(1:length(list_rda), function(i) {
        study_parts <- unlist(strsplit(study_names[i], "__"))
        res <- geneSurvCont(
            dat.icb = list_rda[[i]],
            time.censor = 36,
            missing.perc = 0.5,
            const.int = 0.001,
            n.cutoff = 15,
            feature = "${params.gene_name}",
            study = study_names[i],
            surv.outcome = "OS",
            cancer.type = study_parts[2],
            treatment = study_parts[3]
        )

        # Remove rows with NA coefficients from the results
        res <- res[!is.na(res\$Coef), ]
        res
    })

    # Combine results into one data frame
    cox_os <- do.call(rbind, cox_os)

    # Remove rows without a gene name
    cox_os <- cox_os[!is.na(cox_os\$Gene), ]

    # Adjust p-values for multiple comparisons using the Benjamini-Hochberg method
    cox_os\$FDR <- p.adjust(cox_os\$Pval, method = "BH")

    # Save the results to a CSV file
    write.csv(cox_os, file = "${params.study_id}_meta_analysis_os.csv", row.names = FALSE)
    """
}

//Aggregating Associations (response, R/NR) through Meta-analysis (Pan-cancer)

workflow {

    // Load the specified RDA data file
    icb_dat = file(params.rda_file)

    ////////// Meta Analysis Section //////////

    // Extract expression and clinical data to CSV files
    extracted_data = LoadAndExtractData(icb_dat)

    // Perform gene association analysis
    GeneAssociationOS(extracted_data[0], extracted_data[1]) // expr.csv and clin.csv
    GeneAssociationPFS(extracted_data[0], extracted_data[1])
    LogisticRegression(extracted_data[0], extracted_data[1])


     ////////// Meta Analysis Section //////////
    
    // Load RDA files and extract data
    all_rdas = LoadAllData(params.rda_files_dir)

    // Perform Meta Analysis(OS) using loaded RDA data
    MetaAnalysisOS(all_rdas)
}




