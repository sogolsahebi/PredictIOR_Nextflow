#!/bin/bash -ue
Rscript -e '
library(SummarizedExperiment)
library(dplyr)

# Load the .rda file and assume it directly loads into 'se'
load("ICB_Wolf__Breast__IO+chemo.rda")
se <- get(ls()[1])  # Get the first (or the only) object loaded, assumed to be SummarizedExperiment

# Extract expression data and clinical data
expr <- assay(se)
clin <- as.data.frame(colData(se))

# Source your gene survival analysis script
source("R/geneSurvCont.R")

# Apply gene survival analysis
cox_os <- lapply(1:100, function(i) {
  res <- geneSurvCont(dat.icb = expr,
                      clin = clin,
                      time.censor = 36,
                      missing.perc = 0.5,
                      const.int = 0.001,
                      n.cutoff = 15,
                      feature = rownames(expr)[i],
                      study = "ICB_Ravi__Lung__PD-(L)1",
                      surv.outcome = "OS",
                      cancer.type = "Lung",
                      treatment = "PD-(L)1")

  if (!is.null(res) && !is.na(res$Coef)) {
    return(res)
  } else {
    return(NULL)
  }
})

# Filter NULL results and combine
cox_os <- do.call(rbind, cox_os[!sapply(cox_os, is.null)])
cox_os$FDR <- p.adjust(cox_os$Pval, method = "BH")

# Save results to CSV file
write.csv(cox_os[order(cox_os$FDR, decreasing = FALSE), ], "results.csv", row.names = FALSE)
'
