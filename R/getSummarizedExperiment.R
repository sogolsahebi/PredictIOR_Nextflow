# Goal: 
#- Function to convert a MAE to multiple SummarizedExperiments based on all unique treatment types

# Load required packages
library(SummarizedExperiment)

getSummarizedExperiment <- function(maepath, outputPath) {
  mae <- readRDS(maepath)
  
  #TODO: iprove this: dat_type_ix <- which(dat_type %in% names(dat_icb))
  expr <- assays(mae)[["expr"]]
  clin <- data.frame(colData(mae))
  annot <- data.frame(rowData(mae@ExperimentList$expr))
  
  annot <- annot[annot$gene_type == "protein_coding", , drop=FALSE]
  remove_Par_Y <- grep("PAR_Y", rownames(annot))
  if(length(remove_Par_Y) > 0){
    annot <- annot[-remove_Par_Y, ]
  }
  
  annot <- annot[order(rownames(annot)), , drop=FALSE]
  annot <- annot[!duplicated(annot$gene_name), , drop=FALSE]
  expr <- expr[rownames(annot), , drop=FALSE]
  rownames(expr) <- rownames(annot) <- annot$gene_name
  
  uniqueTreatments <- unique(clin$treatment) # e.g.: "PD-1/PD-L1", "IO+chemo", "IO+combo"
  seList <- list()
  
  for (treatment in uniqueTreatments) {
    clinSubset <- clin[clin$treatment == treatment, ]
    exprSubset <- expr[, rownames(clinSubset)]
    
    print(paste("Processing treatment:", treatment))
    
    se <- SummarizedExperiment(assays = list(gene_expression = exprSubset),
                               rowData = annot, colData = clinSubset)
    
    # Adjusted fileName creation
    fileName <- sprintf("%s/ICB_Ravi__%s__%s.rda", outputPath, clinSubset$cancer_type[1], gsub("[+/]", "_", treatment))
    
    save(se, file=fileName)
    seList[[treatment]] <- se
  }
  
  return(seList)
}

outputPath <- "results/"
maepath <- "data/ICB_Ravi.rds"
mae_SEs <- getSummarizedExperiment(maepath, outputPath)

