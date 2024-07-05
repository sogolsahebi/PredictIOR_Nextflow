#!/usr/bin/env nextflow
nextflow.enable.dsl=2

// Workflow configuration
params.study_id = 'ICB_small_Mariathasan'
params.data_dir = './ICB_data'
params.sig_data_dir = './SIG_data'
params.out_dir = './output'

cancer_type = 'Bladder' // Set cancer type name for the dataset
treatment = 'PD-1/PD-L1' // Set treatment type name for the dataset

/*
Load Immunotherapy RNA Signatures
We evaluate the reproducibility of a compendium of more than 50
RNA signatures described as immunotherapy biomarkers. The following 
published immunotherapy RNA signatures are available on the GitHub 
data repository. To compute the signatures’ scores, different methods 
are used, including weighted mean, single sample gene set enrichment 
analysis (ssGSEA), Gene Set Variation Analysis (GSVA) PMID 23323831, 
and specific computational algorithms.

links:
- GitHub : https://github.com/bhklab/SignatureSets/tree/main
*/

/*
========================================================
SECTION: Signature Score Computation + Assosiation
========================================================

RNA expression signatures:
- Genes without assigned weights are computed using Gene Set Variation Analysis (GSVA) enrichment score with the R package [REF].
- Signatures with specific weights (+1 for up-regulated, -1 for down-regulated genes) use the weighted mean expression approach.
- For RNA immunotherapy signatures requiring specific computational algorithms (PredictIO), original publication methods are employed.

For each dataset:
- RNA signatures are computed if at least 80% of their genes are present.
- Z-score transformation is applied to genes of each RNA signature, before GSVA or weighted mean computation, and after specific immunotherapy signature computation.
*/

/*
---------------------------------------------------------
SUBSECTION: GSVA signature score computation - GSVA approach
---------------------------------------------------------

(PMID 25594174) is a RNA signature capturing the local immune cytolytic
(CYT) activity. To get the signature score, the GSVA approach is applied.
*/

params.sig_name_GSVA = 'CYT_Rooney'

process GeneSigScore_GSVA {
    tag "${params.study_id}_GSVA"
    container 'nextflow-env:latest'
    publishDir "${params.out_dir}", mode: 'copy'

    input:
    path icb_rda_file
    path sig_rda_file

    output:
    path "${params.study_id}_GeneSigScore_GSVA.csv"

    script:
    """
    #!/usr/bin/env Rscript
    library(SummarizedExperiment)
    library(PredictioR)
    library(GSVA)

    # Load icb rda file
    load("${icb_rda_file}")
    expr <- assay(${params.study_id})
    clin <- as.data.frame(colData(${params.study_id}))

    # Load signature genes
    load("${sig_rda_file}") #this gives us 'sig'

    geneSigScore <- geneSigGSVA(dat.icb = expr,
                                clin = clin,
                                sig = sig,
                                sig.name = "${params.sig_name_GSVA}",
                                missing.perc = 0.5,
                                const.int = 0.001,
                                n.cutoff = 15,
                                sig.perc = 0.8, 
                                study = "${params.study_id}")

    geneSigScore <- as.data.frame(geneSigScore)
    # Save as CSV file
    write.csv(geneSigScore, file = "${params.study_id}_GeneSigScore_GSVA.csv", row.names = FALSE)
    """
}

/*
---------------------------------------------------------
SUBSECTION: Weighted mean signature score computation
---------------------------------------------------------
*/

params.sig_name_WeightedMean = 'EMT_Thompson'

process GeneSigScore_WeightedMean {
    tag "${params.study_id}_WeightedMean"
    container 'nextflow-env:latest'
    publishDir "${params.out_dir}", mode: 'copy'

    input:
    path icb_rda_file
    path sig_rda_file

    output:
    path "${params.study_id}_GeneSigScore_WeightedMean.csv"

    script:
    """
    #!/usr/bin/env Rscript
    library(SummarizedExperiment)
    library(PredictioR)

    # Load icb rda file
    load("${icb_rda_file}")
    expr <- assay(${params.study_id})
    clin <- as.data.frame(colData(${params.study_id}))

    # Load signature genes
    load("${sig_rda_file}") # this gives us 'sig'

    geneSigScore <- geneSigMean(dat.icb = expr,
                                        clin = clin,
                                        sig = sig,
                                        sig.name = "${params.sig_name_WeightedMean}",
                                        missing.perc = 0.5,
                                        const.int = 0.001,
                                        n.cutoff = 15,
                                        sig.perc = 0.8, 
                                        study = "${params.study_id}")

    geneSigScore <- as.data.frame(geneSigScore)
    # Save as CSV file
    write.csv(geneSigScore, file = "${params.study_id}_GeneSigScore_WeightedMean.csv", row.names = FALSE)
    """
}

/*
----------------------------------------------------------------------
SUBSECTION: Specific algorithm (PredictIO) signature score computation
----------------------------------------------------------------------
*/

//As an example, ADO (PMID 31953314) is a RNA signature capturing the adenosine (ADO) 
// pathway activity. This adenosine signature was significantly associated with reduced 
// efficacy of anti-PD1 therapy in published cohorts. To get the signature score, the GSVA 
//approach is applied.
//(PMID 36055464) is a RNA signature that demonstrated better and more consistent ability to predict immunotherapy response as compared to other signatures.
// use an rds file Signature_data are from Github

//- link: https://github.com/bhklab/SignatureSets/tree/main/data.

//Signature_data = file("./Bareche/PredictIO_Bareche.rda")
params.sig_name_PredictIO = 'PredictIO_Bareche'

process GeneSigScore_PredictIO {
    tag "${params.study_id}_PredictIO"
    container 'nextflow-env:latest'
    publishDir "${params.out_dir}", mode: 'copy'

    input:
    path icb_rda_file
    path sig_rda_file

    output:
    path "${params.study_id}_GeneSigScore_PredictIO.csv"

    script:
    """
    #!/usr/bin/env Rscript
    library(SummarizedExperiment)
    library(PredictioR)
    library(GSVA)

    # Load icb rda file
    load("${icb_rda_file}")
    expr <- assay(${params.study_id})
    clin <- as.data.frame(colData(${params.study_id}))

    # Load signature genes
    load("${sig_rda_file}") # this gives us 'sig'

    geneSigScore_PredictIO <- geneSigPredictIO(dat.icb = expr,
                                               clin = clin, 
                                               sig = sig,
                                               sig.name = "${params.sig_name_PredictIO}",
                                               missing.perc = 0.5,
                                               const.int = 0.001,
                                               n.cutoff = 15,
                                               sig.perc = 0.8, 
                                               study = "${params.study_id}")

    # Convert to dataframe
    geneSigScore_PredictIO <- as.data.frame(t(geneSigScore_PredictIO))

    # Save as CSV file
    write.csv(geneSigScore_PredictIO, file = "${params.study_id}_GeneSigScore_PredictIO.csv", row.names = FALSE)
    """
}


/*
-----------------------------------------------------------------------------
SUBSECTION: Association of RNA Signatures with Immunotherapy Response &
            Aggregating Associations(OS) through Meta-analysis (Pan-cancer)
-----------------------------------------------------------------------------

For signature, we assess the association of signature with immunotherapy OS.The associations across datasets are aggregated using the meta-analysis approach

all icb rda files are comming from:
- GitHub: https://github.com/bhklab/PredictioR/tree/main/data ,Also loacted in './ICB_data'  directory
- set theis cancer_types and treatments types repesctively
*/

// Set the cancer type and treatment for each dataset respectively.
params.cancer_types = '["Melanoma", "Pancreas", "Melanoma", "Bladder", "Melanoma", "Melanoma", "Kidney", "Melanoma"]'
params.treatment_types = '["PD-1/PD-L1", "PD-1/PD-L1", "PD-1/PD-L1", "PD-1/PD-L1", "CTLA4", "IO+combo", "PD-1/PD-L1", "CTLA4"]'
params.gene_name = "CXCL9"

process RNAGeneSigAssociation {
    container 'nextflow-env:latest'
    publishDir "${params.out_dir}", mode: 'copy'

    input:
    path icb_rda_file
    path sig_rda_file

    output:
    path "${params.study_id}_GeneSigAssociationOS.csv"
    path "${params.study_id}_signature_Associations(OS)_pancancer.csv"

    script:
    """
    #!/usr/bin/env Rscript
    library(SummarizedExperiment)
    library(PredictioR)
    library(survival)
    library(meta)
    library(GSVA)
    library(jsonlite) 

    load("${sig_rda_file}") # this gives us 'sig'

    # Load all .rda files 
    lapply(list.files(path = '${icb_rda_file}', pattern = '*.rda', full.names = TRUE), function(file) {
        load(file, envir = .GlobalEnv)
    })

    # Create a list of the loaded objects using their actual names
    loaded_objects <- ls(pattern = "^ICB_small_")
    expr <- mget(loaded_objects, envir = .GlobalEnv)

    # Define the cancer types and treatment types vectors
    cancer_type <- fromJSON('${params.cancer_types}')
    treatment_type <- fromJSON('${params.treatment_types}')

    geneSig_os <- lapply(1:length(expr), function(k){
       
       geneSig <- geneSigGSVA(dat.icb = expr[[k]], 
                              sig = sig,
                              sig.name = "${params.sig_name_GSVA}",
                              missing.perc = 0.5,
                              const.int = 0.001, 
                              n.cutoff = 15,
                              sig.perc = 0.8, 
                              study = names(expr)[k])
       
       if(sum(!is.na(geneSig)) > 0){
         
         res <- geneSigSurvDicho(dat.icb = expr[[k]],
                                 geneSig = geneSig[1,],
                                 time.censor = 24,
                                 n.cutoff = 15,
                                 study =  names(expr)[k],
                                 surv.outcome = "OS",
                                 n0.cutoff = 5,
                                 n1.cutoff = 5,
                                 sig.name = "${params.sig_name_GSVA}",
                                 method = 'median',
                                 var.type = FALSE,
                                 cancer.type = cancer_type[k],
                                 treatment = treatment_type[k])
         
       } else {
         
         res <- data.frame( Outcome = "OS",
                            Gene = NA, 
                            Study = names(expr)[k],
                            Coef = NA,
                            SE = NA,
                            N = NA,
                            Pval = NA,
                            Cancer_type= NA,
                            Treatment = NA) 
       }
       
       rownames(res) <- NULL
       
       res
       
     })
     
    geneSig_os <- do.call(rbind, geneSig_os)
    geneSig_os\$FDR <- p.adjust(geneSig_os\$Pval, method = "BH")
    geneSig_os <- geneSig_os[order(geneSig_os\$Pval, decreasing = FALSE), ]


    # Now Aggregating Associations(OS) through Meta-analysis (Pan-cancer)

    res <- metafun(coef = geneSig_os\$Coef, 
               se = geneSig_os\$SE, 
               study  = geneSig_os\$Study, 
               pval = geneSig_os\$Pval, 
               n = geneSig_os\$N, 
               cancer.type = geneSig_os\$Cancer_type,
               treatment = geneSig_os\$Treatment,
               feature = unique(geneSig_os\$Gene),
               cancer.spec = FALSE,
               treatment.spec = FALSE)

    write.csv(geneSig_os, file = "${params.study_id}_GeneSigAssociationOS.csv", row.names = FALSE)
    write.csv(res, file = "${params.study_id}_signature_Associations(OS)_pancancer.csv", row.names = FALSE)
    """
}

workflow {
    /*
    ========================================================
    SECTION: Signature Score Computation
    ========================================================
    */
    // Specify input files
    icb_rda_file = file("${params.data_dir}/${params.study_id}.rda")
    sig_rda_file_GSVA = file("${params.sig_data_dir}/${params.sig_name_GSVA}.rda")
    sig_rda_file_WeightedMean = file("${params.sig_data_dir}/${params.sig_name_WeightedMean}.rda")
    sig_rda_file_PredictIO = file("${params.sig_data_dir}/${params.sig_name_PredictIO}.rda")
    sig_data_dir = file(params.sig_data_dir)
    icb_all_rdas = file(params.data_dir)

    // Call the GSVA process
    GeneSigScore_GSVA(icb_rda_file, sig_rda_file_GSVA)
    
    // Call the weighted mean process
    GeneSigScore_WeightedMean(icb_rda_file, sig_rda_file_WeightedMean)
    
    // Call the PredictIO process
    GeneSigScore_PredictIO(icb_rda_file, sig_rda_file_PredictIO)

    /*
    =========================================================================
    SECTION: Signature Association: 
             Association of RNA Signatures with Immunotherapy Response &
             Aggregating Associations(OS) through Meta-analysis (Pan-cancer)
    ===========================================================================
    */

    
    RNAGeneSigAssociation(icb_all_rdas , sig_rda_file_GSVA)
}
