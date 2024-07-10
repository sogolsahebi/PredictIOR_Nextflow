#!/usr/bin/env nextflow
nextflow.enable.dsl=2

// Workflow configuration
params.study_id = 'ICB_small_Mariathasan' // Study identifier, should match with .rda filenames
params.icb_data_dir = './ICB_data' // Directory containing ICB data files
params.sig_data_dir = './SIG_data' // Directory containing signature data files
params.sig_summery_dir = './sig_summery_info' // Directory containing signature information summary 
params.out_dir = './output' // Output directory for results
params.cancer_type = 'Bladder' // Type of cancer being studied
params.treatment = 'PD-1/PD-L1' // Type of treatment applied in the study

log.info """
P R E D I C T I O - N F   P I P E L I N E (Signature Level Analysis)
====================================================================
Study ID             : ${params.study_id}
ICB Data Directory   : ${params.icb_data_dir}
SIG Data Directory   : ${params.sig_data_dir}
SIG info summary     : ${params.sig_summery_dir}
Output Directory     : ${params.out_dir}
Cancer Type          : ${params.cancer_type}
Treatment            : ${params.treatment}
""".stripIndent()

/*

We evaluate the reproducibility of a compendium of more than 50
RNA signatures described as immunotherapy biomarkers. The following 
published immunotherapy RNA signatures are available on the GitHub 
data repository. To compute the signatures’ scores, different methods 
are used, including weighted mean, single sample gene set enrichment 
analysis (ssGSEA), Gene Set Variation Analysis (GSVA) PMID 23323831, 
and specific computational algorithms.

links:
- GitHub : https://github.com/bhklab/SignatureSets/tree/main
- RNA signatures: https://github.com/bhklab/SignatureSets/tree/main/data
- signature_information.csv :(Includes information about each signature): https://github.com/bhklab/SignatureSets/blob/main/data-raw/signature_information.csv 

it include signature (name), DNA/RNA (type), association (resistance/sensitivity),
RNA type (e.g., RNA-seq, log CPM), method (e.g., GSVA), cancer type (types of cancer), 
score function (scoring method), PMID (reference ID), curated information (curated: YES/NO), 
immunotherapy (related: Yes/No), GitHub status (Shared/Not)
*/

/*
========================================================
SECTION: Signature Score Computation 
========================================================

RNA expression signatures:
- Genes without assigned weights are computed using Gene Set Variation Analysis (GSVA) enrichment score with the R package [REF].
- Signatures with specific weights (+1 for up-regulated, -1 for down-regulated genes) use the weighted mean expression approach.
- For RNA immunotherapy signatures requiring specific computational algorithms (PredictIO), original publication methods are employed.

For each dataset:
- RNA signatures are computed if at least 80% of their genes are present.
- Z-score transformation is applied to genes of each RNA signature, before GSVA or weighted mean computation, and after specific immunotherapy signature computation.

Methods that are used in this process: GSVA, ssGSEA, Specific Algorithm, Weighted Mean 
*/

process GeneSigScore {
    tag "${params.study_id}"
    container 'bhklab/nextflow-env'
    publishDir "${params.out_dir}", mode: 'copy'

    input:
    path signature_information
    path sig_rda_path
    path icb_rda_path

    output:
    path "${params.study_id}_GeneSigScore.csv"

    script:
    """
    #!/usr/bin/env Rscript
    source('/R/load_libraries.R')

    signature <- read.csv("${signature_information}")
    signature\$Signature <- as.character(signature\$signature)
    signature\$method <- as.character(signature\$method)

    dir_GeneSig <- '${params.sig_data_dir}'
    GeneSig_list <- list.files(path = '${params.sig_data_dir}', pattern = '*.rda', full.names = TRUE)

    load("${icb_rda_path}")

    geneSig.score <- lapply(1:length(GeneSig_list), function(i) {
        load(GeneSig_list[i])
        sig_name <- substr(basename(GeneSig_list[i]), 1, nchar(basename(GeneSig_list[i])) - 4)


        method <- signature[signature\$Signature == sig_name, "method"]

        if (signature[signature\$Signature == sig_name, "method"] == "GSVA") {
            geneSig <- geneSigGSVA(dat.icb = ${params.study_id}, sig = sig, sig.name = sig_name, missing.perc = 0.5, const.int = 0.001, n.cutoff = 15, sig.perc = 0.8, study = "${params.study_id}")
            if (sum(!is.na(geneSig)) > 0) {
                geneSig <- geneSig[1,]
            }
        } else if (signature[signature\$Signature == sig_name, "method"] == "Weighted Mean") {
            geneSig <- geneSigMean(dat.icb = ${params.study_id}, sig = sig, sig.name = sig_name, missing.perc = 0.5, const.int = 0.001, n.cutoff = 15, sig.perc = 0.8, study = "${params.study_id}")
        } else if (signature[signature\$Signature == sig_name, "method"] == "ssGSEA") {
            geneSig <- geneSigssGSEA(dat.icb = ${params.study_id}, sig = sig, sig.name = sig_name, missing.perc = 0.5, const.int = 0.001, n.cutoff = 15, sig.perc = 0.8, study = "${params.study_id}")
            if (sum(!is.na(geneSig)) > 0) {
                geneSig <- geneSig[1,]
            }
        } else if (signature[signature\$Signature == sig_name, "method"] == "Specific Algorithm" & sig_name == "COX-IS_Bonavita") {
            geneSig <- geneSigCOX_IS(dat.icb = ${params.study_id}, sig = sig, sig.name = signature\$Signature[i], missing.perc = 0.5, const.int = 0.001, n.cutoff = 15, sig.perc = 0.8, study = "${params.study_id}")
        } else if (signature[signature\$Signature == sig_name, "method"] == "Specific Algorithm" & sig_name == "IPS_Charoentong") {
            geneSig <- geneSigIPS(dat.icb = ${params.study_id}, sig = sig, sig.name = signature\$Signature[i], missing.perc = 0.5, const.int = 0.001, n.cutoff = 15, sig.perc = 0.8, study = "${params.study_id}")
        } else if (signature[signature\$Signature == sig_name, "method"] == "Specific Algorithm" & sig_name == "PredictIO_Bareche") {
            geneSig <- geneSigPredictIO(dat.icb = ${params.study_id}, sig = sig, sig.name = signature\$Signature[i], missing.perc = 0.5, const.int = 0.001, n.cutoff = 15, sig.perc = 0.8, study = "${params.study_id}")
        } else if (signature[signature\$Signature == sig_name, "method"] == "Specific Algorithm" & sig_name == "IPRES_Hugo") {
            geneSig <- geneSigIPRES(dat.icb = ${params.study_id}, sig = sig, sig.name = signature\$Signature[i], missing.perc = 0.5, const.int = 0.001, n.cutoff = 15, sig.perc = 0.8, study = "${params.study_id}")
        } else if (signature[signature\$Signature == sig_name, "method"] == "Specific Algorithm" & sig_name == "PassON_Du") {
            geneSig <- geneSigPassON(dat.icb = ${params.study_id}, sig = sig, sig.name = signature\$Signature[i], missing.perc = 0.5, const.int = 0.001, n.cutoff = 15, sig.perc = 0.8, study = "${params.study_id}")
        } else if (signature[signature\$Signature == sig_name, "method"] == "Specific Algorithm" & sig_name == "IPSOV_Shen") {
            geneSig <- geneSigIPSOV(dat.icb = ${params.study_id}, sig = sig, sig.name = signature\$Signature[i], missing.perc = 0.5, const.int = 0.001, n.cutoff = 15, sig.perc = 0.8, study = "${params.study_id}")
        }

        if (sum(!is.na(geneSig)) > 0) {
            geneSig <- geneSig
        } else {
            geneSig <- rep(NA, ncol(${params.study_id}))
        }

        geneSig
    })

    geneSig.score <- do.call(rbind, geneSig.score)
    rownames(geneSig.score) <- substr(basename(GeneSig_list), 1, nchar(basename(GeneSig_list)) - 4)

    remove <- which(is.na(rowSums(geneSig.score)))
    if (length(remove) > 0) {
        geneSig.score <- geneSig.score[-remove, ]
    }
    write.csv(geneSig.score, file = "${params.study_id}_GeneSigScore.csv", row.names = FALSE)
    """
}

/*
========================================================
SECTION: Association for "OS" , "PFS" and "Response"
========================================================
*/

// for OS and PFS
process GeneAssociationSurvival {
    tag "${params.study_id}"
    container 'nextflow-env:latest'
    publishDir "${params.out_dir}", mode: 'copy'

    input:
    path icb_rda_path
    path genescore_path
    val survival_type

    output:
    path("${params.study_id}_${survival_type}_GeneSig_association.csv")

    script:
    """
    #!/usr/bin/env Rscript
    source('/R/load_libraries.R')

    load("${icb_rda_path}")
    
    geneSig.score <- read.csv("${genescore_path}")

    res.all <- lapply(1:nrow(geneSig.score), function(k) {
    sig_name <- rownames(geneSig.score)[k]
    
    res <- geneSigSurvCont(
        dat.icb = ${params.study_id},
        geneSig = as.numeric(geneSig.score[k, ]),  
        time.censor = 36,
        n.cutoff = 15,
        study = "${params.study_id}",
        surv.outcome = "${survival_type}",
        sig.name = sig_name,
        cancer.type = "${params.cancer_type}",
        treatment ="${params.treatment}"
    )

    res
    })

    res.all <- do.call(rbind, res.all)
    res.all\$FDR <- p.adjust(res.all\$Pval, method="BH")
    res.all <- res.all[order(res.all\$FDR), ]

    write.csv(res.all, file = "${params.study_id}_${survival_type}_GeneSig_association.csv", row.names = FALSE)
    """
}

// for "Response"(R vs NR)
process GeneAssociationResponse {
    tag "${params.study_id}"
    container 'nextflow-env:latest'
    publishDir "${params.out_dir}", mode: 'copy'

    input:
    path icb_rda_path
    path genescore_path

    output:
    path("${params.study_id}_GeneSig_Response.csv")

    script:
    """
    #!/usr/bin/env Rscript
    source('/R/load_libraries.R')

    load("${icb_rda_path}")
    geneSig.score <- read.csv("${genescore_path}")

    res.logreg <- lapply(1:nrow(geneSig.score), function(k){
    sig_name <- rownames(geneSig.score)[k]
    res <- geneSigLogReg(dat.icb = ${params.study_id},
                        geneSig = as.numeric(geneSig.score[k,]),
                        n.cutoff = 10,
                        study =  "${params.study_id}",
                        sig.name = sig_name,
                        n0.cutoff = 3, 
                        n1.cutoff = 3,
                        cancer.type = "${params.cancer_type}",
                        treatment = "${params.treatment}")
    
    res
    })

    res.logreg <- do.call(rbind, res.logreg)
    res.logreg\$FDR <- p.adjust(res.logreg\$Pval, method="BH")
    # Save as CSV file
    write.csv(res.logreg, file = "${params.study_id}_GeneSig_Response.csv", row.names = FALSE)
    """
}


workflow {


/*
========================================================
 Signature Score Computation 
========================================================
*/ 
    // Specify input files
    icb_rda_path = file("${params.icb_data_dir}/${params.study_id}.rda")
    all_sigs = file(params.sig_data_dir)
    sigs_info_path = file("${params.sig_summery_dir}/signature_information.csv")

gene_sigscore = GeneSigScore(sigs_info_path, all_sigs, icb_rda_path)

/*
========================================================
Sig Association with OS or PFS
========================================================
*/ 

// survival type can be chose a "OS" or "PFS"
survival_type = "OS"
gene_sigscore = gene_sigscore[0]
GeneAssociationSurvival( icb_rda_path, gene_sigscore,survival_type)

/*
========================================================
Sig Association with Response
========================================================
*/ 
GeneAssociationResponse( icb_rda_path, gene_sigscore)

}
