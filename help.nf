/*MetaAnalysis_PerCancer(
        cox_files: all_rda_files.collect { GeneAssociationSurvivalOS(it, clin_file, genes) },
        params.gene_name
    )

    MetaAnalysis_PerTreatment(
        cox_files: all_rda_files.collect { GeneAssociationSurvivalOS(it, clin_file, genes) },
        params.gene_name
    )

    // Perform meta-analysis for PFS
    MetaAnalysis_PanCancer(
        cox_files: all_rda_files.collect { GeneAssociationSurvivalPFS(it, clin_file, genes) },
        params.gene_name
    )

    MetaAnalysis_PerCancer(
        cox_files: all_rda_files.collect { GeneAssociationSurvivalPFS(it, clin_file, genes) },
        params.gene_name
    )

    MetaAnalysis_PerTreatment(
        cox_files: all_rda_files.collect { GeneAssociationSurvivalPFS(it, clin_file, genes) },
        params.gene_name
    )

    // Perform meta-analysis for Response
    MetaAnalysis_PanCancer(
        cox_files: all_rda_files.collect { GeneAssociationResponse(it, clin_file, genes) },
        params.gene_name
    )

    MetaAnalysis_PerCancer(
        cox_files: all_rda_files.collect { GeneAssociationResponse(it, clin_file, genes) },
        params.gene_name
    )

    MetaAnalysis_PerTreatment(
        cox_files: all_rda_files.collect { GeneAssociationResponse(it, clin_file, genes) },
        params.gene_name
    )
    */