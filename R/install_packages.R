#!/usr/bin/env Rscript

# Set the default library path for package installations
.libPaths("~/R/x86_64-pc-linux-gnu-library/4.3")

# List of required packages
packages <- c('survcomp', 'GSVA', 'dplyr', 'meta', 'metafor', 'MultiAssayExperiment', 
              'data.table', 'forestplot', 'ggplot2', 'ggrepel', 'gridExtra', 
              'kableExtra', 'summarytools', 'languageserver')

# Function to check and install missing packages
install_if_missing <- function(packages) {
  for (pkg in packages) {
    if (!require(pkg, character.only = TRUE, quietly = TRUE)) {
      install.packages(pkg, dependencies = TRUE, 
                       repos = 'https://cloud.r-project.org/', 
                       lib = "~/R/x86_64-pc-linux-gnu-library/4.3")
      if (!require(pkg, character.only = TRUE, quietly = TRUE)) {
        cat(paste0("Failed to install ", pkg, "\n"))
      }
    } else {
      cat(paste0(pkg, " is already installed.\n"))
    }
  }
}

# Install missing packages
install_if_missing(packages)
