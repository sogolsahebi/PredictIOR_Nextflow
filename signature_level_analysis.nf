#!/usr/bin/env nextflow
nextflow.enable.dsl=2

// Workflow configuration
params.study_id = 'ICB_small_Padron' // Identifier for the study (name of the study in ./ICB_data directory)
params.icb_data_dir = './ICB_data' // Directory for ICB data files
params.sig_data_dir = './SIG_data' // Directory for signature data files
params.sig_summery_dir = './sig_summery_info' // Directory for signature summary information
params.out_dir = './output/signatures_level_output' // Output directory for results
params.cancer_type = 'Pancreas' // Type of cancer being studied
params.treatment = 'PD-1/PD-L1' // Type of treatment being studied

 /*
Note
Another input example would be:

params.study_id = 'ICB_small_Liu' 
params.cancer_type = 'Melanoma' 
params.treatment = 'PD-1/PD-L1'
/* 

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
Evaluation of over 50 RNA signatures described as immunotherapy biomarkers. 
We use several methods for computing the scores of these signatures:
- GSVA
- ssGSEA
- Weighted mean expression
- Specific computational algorithms (PredictIO)

Reference Links:
- GitHub Repository: https://github.com/bhklab/SignatureSets/tree/main
- RNA signatures: https://github.com/bhklab/SignatureSets/tree/main/data
- Signature information CSV: https://github.com/bhklab/SignatureSets/tree/main/data-raw
*/

/*
SECTION: Signature Score Computation 
- Genes without weights are scored using GSVA enrichment.
- Signatures with weights (+1 for up, -1 for down) use a weighted mean expression approach.
- Specific computational algorithms are used as described in original publications.
- Each dataset undergoes a Z-score transformation post-computation.
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
    write.csv(geneSig.score, file = "${params.study_id}_GeneSigScore.csv", row.names = TRUE)
    """
 
}


/*
=================================
SECTION: Sig Association for "OS" 
=================================
*/


process GeneSig_AssociationOS {
    tag "${params.study_id}"
    container 'nextflow-env:latest'
    publishDir "${params.out_dir}", mode: 'copy'

    input:
    path icb_rda_path
    path genescore_path

    output:
    path("${params.study_id}_os_GeneSig_association.csv")

    script:
    """
    #!/usr/bin/env Rscript
    source('/R/load_libraries.R')

    load("${icb_rda_path}")
    
    geneSig.score <- read.csv("${genescore_path}", row.names = 1)

    res.all <- lapply(1:nrow(geneSig.score), function(k) {
        sig_name <- rownames(geneSig.score)[k]
        geneSig_vector <- as.numeric(geneSig.score[k, ])
        geneSig_vector <- geneSig_vector[!is.na(geneSig_vector)]

        res <- geneSigSurvCont(
            dat.icb = ${params.study_id},
            geneSig = geneSig_vector,
            time.censor = 36,
            n.cutoff = 15,
            study = "${params.study_id}",
            surv.outcome = "OS",
            sig.name = sig_name,
            cancer.type = "${params.cancer_type}",
            treatment = "${params.treatment}"
        )
        
        res
    })

    res.all <- do.call(rbind, res.all)
    res.all\$FDR <- p.adjust(res.all\$Pval, method="BH")
    res.all <- res.all[order(res.all\$FDR), ]

    write.csv(res.all, file = "${params.study_id}_os_GeneSig_association.csv", row.names = TRUE)
    """
}


/*
=============================
SECTION: Association for "PFS" 
=============================
*/

process GeneSig_AssociationPFS {
    tag "${params.study_id}"
    container 'nextflow-env:latest'
    publishDir "${params.out_dir}", mode: 'copy'

    input:
    path icb_rda_path
    path genescore_path

    output:
    path("${params.study_id}_pfs_GeneSig_association.csv")

    script:
    """
    #!/usr/bin/env Rscript
    source('/R/load_libraries.R')

    load("${icb_rda_path}")
    
    geneSig.score <- read.csv("${genescore_path}", row.names = 1)
    

    res.all <- lapply(1:nrow(geneSig.score), function(k) {
    sig_name <- rownames(geneSig.score)[k]
    geneSig_vector <- as.numeric(geneSig.score[k, ])
    geneSig_vector <- geneSig_vector[!is.na(geneSig_vector)]
    
    res <- geneSigSurvCont(
        dat.icb = ${params.study_id},
        geneSig = geneSig_vector,  
        time.censor = 24,
        n.cutoff = 15,
        study = "${params.study_id}",
        surv.outcome = "PFS",
        sig.name = sig_name,
        cancer.type = "${params.cancer_type}",
        treatment ="${params.treatment}"
    )

    res
    })

    res.all <- do.call(rbind, res.all)
    res.all\$FDR <- p.adjust(res.all\$Pval, method="BH")
    res.all <- res.all[order(res.all\$FDR), ]

    write.csv(res.all, file = "${params.study_id}_pfs_GeneSig_association.csv", row.names = TRUE)
    """
}

/*
========================================================
SECTION: Association for Response(R vs NR)
========================================================
*/

//asososiation of signiture with one study 
process GeneSig_AssociationResponse {
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
    geneSig.score <- read.csv("${genescore_path}", row.names = 1)

    res.logreg <- lapply(1:nrow(geneSig.score), function(k){
        sig_name <- rownames(geneSig.score)[k]
        geneSig_vector <- as.numeric(geneSig.score[k, ])
        geneSig_vector <- geneSig_vector[!is.na(geneSig_vector)]

        res <- geneSigLogReg(dat.icb = ${params.study_id},
                            geneSig = geneSig_vector,
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
    write.csv(res.logreg, file = "${params.study_id}_GeneSig_Response.csv", row.names = TRUE)
    """
}

workflow {

/*
========================================================
 Signature Score Computation 
========================================================
*/ 

// set input files
icb_rda_path = file("${params.icb_data_dir}/${params.study_id}.rda")
all_sigs = file(params.sig_data_dir)
sigs_info_path = file("${params.sig_summery_dir}/signature_information.csv")

gene_sigscore_result = GeneSigScore(sigs_info_path, all_sigs , icb_rda_path)


/*
========================================================
 Signature Score Computation 
========================================================
*/

gene_sigscore = gene_sigscore_result[0] 
GeneSig_AssociationOS(icb_rda_path, gene_sigscore)
GeneSig_AssociationPFS(icb_rda_path, gene_sigscore) 
GeneSig_AssociationResponse(icb_rda_path, gene_sigscore)

}