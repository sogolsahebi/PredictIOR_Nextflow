#!/usr/bin/env Rscript
library(SummarizedExperiment)

expr <- read.csv("ICB_Ravi__Lung__PD-L1_expr.csv", row.names = 1)
clin <- read.csv("ICB_Ravi__Lung__PD-L1_clin.csv")

source('/R_scripts/getGeneAssociation.R')
source('/R_scripts/getHR.R')

results <- lapply(1:100, function(i) {
    geneSurvCont(
        dat.icb = expr,
        clin = clin,
        time.censor = 36,
        missing.perc = 0.5,
        const.int = 0.001,
        n.cutoff = 15,
        feature = rownames(expr)[i],
        study = 'ICB_Ravi__Lung__PD-L1',
        surv.outcome = 'OS',
        cancer.type = 'Lung',
        treatment = 'PD-L1'
    )
})

cox_os <- do.call(rbind, results)

# Write the head of the dataframe to a CSV file
write.csv(head(cox_os), file = "ICB_Ravi__Lung__PD-L1_cox_os_head.csv", row.names = TRUE)

# Write column names to a text file
write(colnames(cox_os), file = "ICB_Ravi__Lung__PD-L1_cox_os_colnames.txt")

# Assuming 'Gene' should be a column, handle if not present
if (!"Gene" %in% colnames(cox_os)) {
    stop("Error: 'Gene' column not found in cox_os dataframe")
}



write.csv(cox_os, file = "ICB_Ravi__Lung__PD-L1_cox_os_results.csv", row.names = FALSE)
