#!/usr/bin/env Rscript
source('/R/load_libraries.R')

# Define the directory where your files are located
sig_level_result_dir <- 'signatures_level_output'

# Determine the pattern based on io_outcome
pattern <- ifelse('OS' == 'OS', '_os', ifelse('OS' == 'PFS', '_pfs', '_Response'))

# List all files that contain the pattern in their filenames
sig_files <- list.files(path = sig_level_result_dir, pattern = pattern, full.names = TRUE)

# Read each file and store the results
res <- lapply(sig_files, function(file) {
    read.csv(file)
})

# Combine the results into a single data frame
res <- do.call(rbind, res)
res <- res[!is.na(res$Coef), ]

# Convert to data frame
df <- res
signature <- unique(df$Gene)

# Perform per-cancer meta-analysis on each gene signature
AllGeneSig_meta <- lapply(1:length(signature), function(j) {
    print(j)
    sub_df <- df[df$Gene == signature[j], ]
    if (nrow(sub_df) >= 3) {
        res <- metaPerCanfun(coef = sub_df$Coef,
                             se = sub_df$SE,
                             study  = sub_df$Study,
                             pval = sub_df$Pval,
                             n = sub_df$N,
                             cancer.type = sub_df$Cancer_type,
                             treatment = sub_df$Treatment,
                             cancer.spec = TRUE,
                             feature = unique(sub_df$Gene))

        percan_res <- lapply(1:length(res), function(i) {
            res[[i]]$meta_summery
        })

        percan_res <- do.call(rbind, percan_res)
    } else {
        percan_res <- data.frame(Cancer_type = "Not Applicable",
                                 Gene = signature[j],
                                 Coef = NA,
                                 SE = NA,
                                 CI_lower = NA,
                                 CI_upper = NA,
                                 Pval = NA,
                                 I2 = NA,
                                 Q_Pval = NA)
    }
    percan_res
})

AllGeneSig_meta <- do.call(rbind, AllGeneSig_meta)
AllGeneSig_meta <- AllGeneSig_meta[!is.na(AllGeneSig_meta$Coef), ]

# FDR adjustment
group <- unique(AllGeneSig_meta$Cancer_type)
AllGeneSig_meta <- lapply(1:length(group), function(k) {
    sub_df <- AllGeneSig_meta[AllGeneSig_meta$Cancer_type == group[k], ]
    sub_df$FDR <- p.adjust(sub_df$Pval, method = "BH")
    sub_df
})

AllGeneSig_meta <- do.call(rbind, AllGeneSig_meta)
AllGeneSig_meta <- AllGeneSig_meta[order(AllGeneSig_meta$FDR), ]

# Save the results to a CSV file
write.csv(AllGeneSig_meta, file = "Meta_analysis_Sig_OS_percancer.csv", row.names = FALSE)
