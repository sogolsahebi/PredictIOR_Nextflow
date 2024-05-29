#!/usr/bin/env Rscript

# Load each RDA file and extract data
list_rda <- lapply(list.files(path = 'data', pattern = '*.rda', full.names = TRUE), function(file) {
    load(file)
    dat_icb
})

# Assign study names to the list elements
study_names <- substr(list.files('data'), 5, nchar(list.files('data')) - 4)
names(list_rda) <- study_names
