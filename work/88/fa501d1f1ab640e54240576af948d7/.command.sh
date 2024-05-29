#!/usr/bin/env Rscript
library(SummarizedExperiment)
library(metafor)

# Load each RDA file and extract data
expr <- lapply(rda_files, function(file) {
    load(file)
    dat_icb
})

    #!/usr/bin/env Rscript

# Load each RDA file and extract data
expr <- lapply(list.files(path = 'data', pattern = '.rda', full.names = TRUE), function(file) {
    load(file)
    dat_icb
