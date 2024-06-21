#!/usr/bin/env nextflow
nextflow.enable.dsl=2

params.study_id = 'ICB_small_Mariathasan'
params.data_dir = './data'
params.out_dir = './output'
params.rda_file = "${params.data_dir}/${params.study_id}.rda"

// Extract cancer type and treatment type of this dataset
//cancer_type = params.study_id.split('__')[1]
//treatment_type = params.study_id.split('__')[2]
//println("for ${params.study_id} study has Cancer Type: ${cancer_type}, Treatment Type: ${treatment_type} ")

/*
========================================================
SECTION: Load Immunotherapy Datasets
========================================================
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
========================================================
SECTION: Gene Association Analysis
========================================================
*/

/*
--------------------------------------------------------
SUBSECTION: Overal survival
--------------------------------------------------------
*/

// Gene Association Analysis(OS) for specific gene:
//  Association of CXCL9 gene with OS for selected dataset Mariathasan Bladder
params.gene_name = "CXCL9"

process GeneAssociationOS_SingleGene {
    tag "${params.gene_name}"
    container 'nextflow-env:latest'
    publishDir "${params.out_dir}", mode: 'copy'

    input:
    path expr_csv
    path clin_csv

    output:
    path "${params.study_id}_cox_os_${params.gene_name}.csv"

    script:
    """
    #!/usr/bin/env Rscript
    library(SummarizedExperiment)
    library(PredictioR)
    library(survival)

    expr <- as.matrix(read.csv("${expr_csv}", row.names = 1))
    clin <- read.csv("${clin_csv}", row.names = 1)

    # Ensure cancer_type is extracted correctly
    cancer_type <- names(table(clin\$cancer_type)[table(clin\$cancer_type) >= 15])
       
    # Perform survival analysis
    cox_os <- survCont(
        status = clin\$event_occurred_os,
        time = clin\$survival_time_os,
        time.censor = 36,
        var = as.numeric(expr["${params.gene_name}", ])

    )
    # convert to dataframe
    cox_os <- as.data.frame(t(cox_os))

    write.csv(cox_os, file = "${params.study_id}_cox_os_${params.gene_name}.csv", row.names = FALSE)
    """
}

// This process performs the association analysis of a dichotomous gene with overall survival response to immunotherapy.
process GeneAssociationOS_DichotomousGene{
    tag "${params.gene_name}"
    container 'nextflow-env:latest'
    publishDir "${params.out_dir}", mode: 'copy'

    input:
    path expr_csv
    path clin_csv

    output:
    path "${params.study_id}_dicho_os_${params.gene_name}.csv"

    script:
    """
    #!/usr/bin/env Rscript
    library(SummarizedExperiment)
    library(PredictioR)
    library(survival)

    expr <- as.matrix(read.csv("${expr_csv}", row.names = 1))
    clin <- read.csv("${clin_csv}", row.names = 1)

    # Ensure cancer_type is extracted correctly
    cancer_type <- names(table(clin\$cancer_type)[table(clin\$cancer_type) >= 15])
       
    # Perform survival analysis
    dicho <- survDicho( status = clin\$event_occurred_os ,
           time = clin\$survival_time_os,
           time.censor= 36,
           var = as.numeric(expr["${params.gene_name}", ]),
           n0.cutoff = 5,
           n1.cutoff = 5,
           method = "median",
           var.type = FALSE)
    
    
    # convert to dataframe
    dicho <- as.data.frame(t(dicho))

    write.csv(dicho, file = "${params.study_id}_dicho_os_${params.gene_name}.csv", row.names = FALSE)
    """
}

//Association of genes in PredictIO Bareche RNA signature 36055464 (https://pubmed.ncbi.nlm.nih.gov/36055464/) with immunotherapy OS
process GeneAssociationOS_genes {
    container 'nextflow-env:latest'
    publishDir "${params.out_dir}", mode: 'copy'

    input:
    path expr_file
    path clin_file
    path bareche_file

    output:
    path "${params.study_id}_cox_os_barechegenes.csv"

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
    clin <- read.csv("${clin_file}")
    cancer_type <- names(table(clin\$cancer_type)[table(clin\$cancer_type) >= 15])

    # load "PredictIO_Bareche.rda"
    load("${bareche_file}")

    genes <- PredictIO_Bareche\$gene_name

    cox_os <- geneSurvCont(
        dat.icb = expr,
        clin = clin,
        time.censor = 36,
        missing.perc = 0.5,
        const.int = 0.001,
        n.cutoff = 15,
        feature = genes,
        study = "${params.study_id}",
        surv.outcome = 'OS',
        cancer.type = cancer_type,
        treatment = 'PD-1/PD-L1'
    )

    # Additional filtering 
    cox_os <- cox_os[order(cox_os\$FDR <- p.adjust(cox_os\$Pval, method = "BH")), ]
    write.csv(cox_os, file = "${params.study_id}_cox_os_barechegenes.csv", row.names = FALSE)
    """
}

/*
--------------------------------------------------------
SUBSECTION: Immunotherapy response (R vs NR).
--------------------------------------------------------
*/

process GeneAssociationResponse_RvsNR {
    tag "${params.study_id}"
    container 'nextflow-env:latest'
    publishDir "${params.out_dir}", mode: 'copy'

    input:
    path expr_file
    path clin_file
    path bareche_file

    output:
    path "${params.study_id}_Response_RvsNR.csv"

    script:
    """
    #!/usr/bin/env Rscript
    library(SummarizedExperiment)
    library(PredictioR)

    expr <- read.csv("${expr_file}", row.names = 1)
    clin <- read.csv("${clin_file}")
    cancer_type <- names(table(clin\$cancer_type)[table(clin\$cancer_type) >= 15])

    # load "PredictIO_Bareche.rda"
    load("${bareche_file}")

    genes <- PredictIO_Bareche\$gene_name
    
    logreg <- geneLogReg(dat.icb = expr,
                 clin = clin,
                     missing.perc = 0.5,
                     const.int = 0.001,
                     n.cutoff = 15,
                     feature = genes,
                     study = "${params.study_id}", 
                     n0.cutoff = 10,
                     n1.cutoff = 10,
                     cancer.type = cancer_type,
                     treatment = 'PD-1/PD-L1')

    # Adjust P-values and sort by FDR
    logreg <- logreg[order(logreg\$FDR <- p.adjust(logreg\$Pval, method = "BH")), ] 

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
------------------------------------------------------------------------------
    A.SUBSECTION: Aggregating Associations through Meta-analysis (Pan-cancer) &
    B.SUBSECTION: Aggregating Associations through Meta-analysis (Per-cancer)
    C.Aggregating Associations through Meta-analysis (Per-treatment)
--------------------------------------------------------------------------------
*/
//set "CXCL9" gene for Meta Analysis, 
params.meta_analysis_gene_name = "CXCL9"
process MetaAnalysisOS{
    container 'nextflow-env:latest'
    publishDir "${params.out_dir}", mode: 'copy'

    input:
    path rda_files_path

    output:
    path "${params.study_id}_meta_analysis_os__${params.meta_analysis_gene_name}_assosiation.csv"
    path "${params.study_id}_meta_analysis_os__${params.meta_analysis_gene_name}_pancancer.csv"
    path "${params.study_id}_meta_analysis_os__${params.meta_analysis_gene_name}_percancer.csv"
    path "${params.study_id}_meta_analysis_os__${params.meta_analysis_gene_name}_pertreatment.csv"

    script:
    """
    #!/usr/bin/env Rscript
    library(SummarizedExperiment)
    library(PredictioR)
    library(survival)
    library(meta)

    list_rda <- lapply(list.files(path = '${rda_files_path}', pattern = '*.rda', full.names = TRUE), function(file) {
        load(file)
    })

    expr <- list('ICB_Liu' = ICB_small_Liu, 'ICB_Padron' = ICB_small_Padron, 'ICB_Hugo' = ICB_small_Hugo, 
              'ICB_Mariathasan' = ICB_small_Mariathasan, 'ICB_Nathanson' = ICB_small_Nathanson, 
              'ICB_Riaz' = ICB_small_Riaz, 'ICB_Miao' = ICB_small_Miao, 'ICB_Van_Allen' = ICB_small_Van_Allen)

    cancer_type <- c('Melanoma', 'Pancreas', 'Melanoma', 'Bladder', 'Melanoma', 'Melanoma', 'Kidney', 'Melanoma')
    treatment_type <- c('PD-1/PD-L1', 'PD-1/PD-L1', 'PD-1/PD-L1', 'PD-1/PD-L1', 'CTLA4', 'IO+combo', 'PD-1/PD-L1', 'CTLA4')

    # Apply a function over the loaded datasets to perform survival analysis
    assoc.res <- lapply(1:length(expr), function(k){
        geneSurvCont(dat.icb = expr[[k]],
                     time.censor = 36,
                     missing.perc = 0.5,
                     const.int = 0.001,
                     n.cutoff = 15,
                     feature = "${params.meta_analysis_gene_name}",
                     study = names(expr)[k],
                     surv.outcome = 'OS',
                     cancer.type = cancer_type[k],
                     treatment = treatment_type[k])
    })

    # Combine results into one data frame
    assoc.res <- do.call(rbind, assoc.res)

    # Adjust p-values for multiple comparisons using the Benjamini-Hochberg method
    assoc.res\$FDR <- p.adjust(assoc.res\$Pval, method = "BH")

    # meta-analysis for a gene across datasets
    res_meta_pancancer <- metafun(coef = assoc.res\$Coef, 
                        se = assoc.res\$SE,
                        study  = assoc.res\$Study, 
                        pval = assoc.res\$Pval, 
                        n = assoc.res\$N,
                        cancer.type = assoc.res\$Cancer_type,
                        treatment = assoc.res\$Treatment,
                        feature = "${params.gene_name}", 
                        cancer.spec = FALSE, 
                        treatment.spec = FALSE)
        
    # meta-analysis results
    res_meta_pancancer <- data.frame(res_meta_pancancer)

    res_meta_percancer <- metaPerCanfun(coef = assoc.res\$Coef, 
                          se = assoc.res\$SE,
                          study  = assoc.res\$Study, 
                          pval = assoc.res\$Pval, 
                          n = assoc.res\$N,
                          cancer.type = assoc.res\$Cancer_type,
                          treatment = assoc.res\$Treatment,
                          feature = "${params.gene_name}", 
                          cancer.spec = TRUE)

   res_meta_percancer  <- data.frame(rbind(res_meta_percancer \$Melanoma\$meta_summery,res_meta_percancer\$Other\$meta_summery))

   # treatment specific meta-analysis for a gene across datasets
   res_meta_pertreatment <- metaPerTreatmentfun(coef = assoc.res\$Coef, 
                                    se = assoc.res\$SE,
                                    study  = assoc.res\$Study, 
                                    pval = assoc.res\$Pval, 
                                    n = assoc.res\$N,
                                    cancer.type = assoc.res\$Cancer_type,
                                    treatment = assoc.res\$Treatment,
                                    feature = "${params.gene_name}", 
                                    treatment.spec = TRUE)
        
    # summary of meta-analysis results for each treatment types
    res_meta_pertreatment <- rbind(res_meta_pertreatment\$`PD-1/PD-L1`\$meta_summery,res_meta_pertreatment\$Other\$meta_summery)

    # Save the results to a CSV file
    write.csv(assoc.res, file = "${params.study_id}_meta_analysis_os__${params.meta_analysis_gene_name}_assosiation.csv", row.names = FALSE)
    write.csv(res_meta_pancancer, file = "${params.study_id}_meta_analysis_os__${params.meta_analysis_gene_name}_pancancer.csv", row.names = FALSE)
    write.csv(res_meta_percancer , file = "${params.study_id}_meta_analysis_os__${params.meta_analysis_gene_name}_percancer.csv", row.names = FALSE)
    write.csv(res_meta_pertreatment , file = "${params.study_id}_meta_analysis_os__${params.meta_analysis_gene_name}_pertreatment.csv", row.names = FALSE)

    """
}

    /*
========================================================
    SECTION: Signature Score Computation
========================================================
    */






workflow {
/*
========================================================
SECTION: Load RDA Data and Extract Expression and Clinical Data
========================================================
*/

    // Set the specified study id and data directory
    icb_dat = file(params.rda_file)

    // Extract expression and clinical data to CSV files
    extracted_data = LoadAndExtractData(icb_dat) // extracted_data[0] -> expr.csv, extracted_data[1] -> clin.csv

    /*
========================================================
    SECTION: Gene Association Analysis
========================================================
    */
    // Perform gene association analysis

    /*
    --------------------------------------------------------
    A.SUBSECTION: Overal survival
    --------------------------------------------------------
    */

    GeneAssociationOS_SingleGene(extracted_data[0], extracted_data[1]) // expr.csv and clin.csv
    GeneAssociationOS_DichotomousGene(extracted_data[0], extracted_data[1])

    /*
    --------------------------------------------------------
    B.SUBSECTION: Immunotherapy response (R vs NR).
    --------------------------------------------------------
    */
    // Replace with Bareche data
    Bareche = file("./Bareche/PredictIO_Bareche.rda")
    GeneAssociationOS_genes(extracted_data[0], extracted_data[1], Bareche )
    //GeneAssociationPFS(extracted_data[0], extracted_data[1])
    GeneAssociationResponse_RvsNR (extracted_data[0], extracted_data[1], Bareche )

    /*

========================================================
    SECTION: Meta Analysis -OS
========================================================
    */

    //params.meta_analysis_gene_name = "CXCL9" // Set a gene name as default
    // Load RDA files and extract data
    all_rdas = file(params.data_dir)

    // Perform Meta Analysis (OS) using loaded RDA data
    MetaAnalysisOS(all_rdas)

    /*
========================================================
    SECTION: Signature Score Computation
========================================================
    */

    // Placeholder for future signature analysis steps
    // Add relevant functions and parameters here
}
