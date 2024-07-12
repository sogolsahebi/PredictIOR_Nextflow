#!/usr/bin/env Rscript
source('/R/load_libraries.R')

# Define the directory where your files are located
sig_level_result_dir <- 'signatures_level_output'

if ('OS' == 'OS'){
# List all files that contain '_os' in their filenames
sig_files <- list.files(path = sig_level_result_dir, pattern = '_os', full.names = TRUE)

# Combine file paths
dir <- file.path(sig_level_result_dir, sig_files)

# Load each file and store the results
res.os <- lapply(1:length(dir), function(k) {
  load(dir[k])
  res.os
})

# Combine the results into a single data frame
res.os <- do.call(rbind, res.os)
res.os <- res.os[!is.na(res.os$Coef), ]

# Convert to data frame
df <- res.os
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

}

# Combine meta-analysis results
AllGeneSig_meta <- do.call(rbind, AllGeneSig_meta)
AllGeneSig_meta <- AllGeneSig_meta[!is.na(AllGeneSig_meta$Coef), ]
AllGeneSig_meta$FDR <- p.adjust(AllGeneSig_meta$Pval, method = "BH")
AllGeneSig_meta <- AllGeneSig_meta[order(AllGeneSig_meta$FDR), ]

# Save the results to a CSV file
write.csv(AllGeneSig_meta, file = "Meta_analysis_Sig_OS_pancancer.csv", row.names = FALSE)
