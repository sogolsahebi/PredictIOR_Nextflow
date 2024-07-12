#!/usr/bin/env Rscript
source('/R/load_libraries.R')

# Define the directory where your files are located
sig_level_result_dir <- 'signatures_level_output'

# Determine the pattern based on io_outcome
pattern <- ifelse('Response' == 'OS', '_os', ifelse('Response' == 'PFS', '_pfs', '_Response'))

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

# Perform meta-analysis on each gene signature
AllGeneSig_meta <- lapply(1:length(signature), function(j) {
    print(j)
    res <- metafun(coef = df[df$Gene == signature[j], "Coef"],
                   se = df[df$Gene == signature[j], "SE"],
                   study  = df[df$Gene == signature[j], "Study"],
                   pval = df[df$Gene == signature[j], "Pval"],
                   n = df[df$Gene == signature[j], "N"],
                   cancer.type = df[df$Gene == signature[j], "Cancer_type"],
                   treatment = df[df$Gene == signature[j], "Treatment"],
                   cancer.spec = FALSE,
                   treatment.spec = FALSE,
                   feature = unique(df[df$Gene == signature[j], "Gene"]))

    res$meta_summery
})

# Combine meta-analysis results
AllGeneSig_meta <- do.call(rbind, AllGeneSig_meta)
AllGeneSig_meta <- AllGeneSig_meta[!is.na(AllGeneSig_meta$Coef), ]
AllGeneSig_meta$FDR <- p.adjust(AllGeneSig_meta$Pval, method = "BH")
AllGeneSig_meta <- AllGeneSig_meta[order(AllGeneSig_meta$FDR), ]

# Save the results to a CSV file
write.csv(AllGeneSig_meta, file = "Meta_analysis_Sig_Response_pancancer.csv", row.names = FALSE)
save(AllGeneSig_meta, file = "meta_pan_Response.RData")