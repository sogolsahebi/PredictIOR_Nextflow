#!/usr/bin/env Rscript
source('/R/load_libraries.R')

# Load all .rda files 
lapply(list.files(path = 'ICB_data', pattern = '*.rda', full.names = TRUE), function(file) { 
    load(file, envir = .GlobalEnv) 
})

# Create a list of the loaded objects using their actual names
loaded_objects <- basename(list.files(path = 'ICB_data', pattern = '*.rda', full.names = TRUE))
loaded_objects <- substr(loaded_objects, 1, nchar(loaded_objects) - 4) # remove .rda from name
expr <- mget(loaded_objects, envir = .GlobalEnv)

# Print structure of each loaded dataset for debugging
print("Loaded datasets:")
for (name in names(expr)) {
    print(paste("Dataset:", name))
    print("Variable names:")
    print(names(expr[[name]]))
    print("First few rows:")
    print(head(expr[[name]]))
}

# Define the cancer types and treatment types vectors
cancer_types <- fromJSON('["Melanoma", "Pancreas", "Melanoma", "Bladder", "Melanoma", "Melanoma", "Kidney", "Melanoma"]')
treatment_types <- fromJSON('["PD-1/PD-L1", "PD-1/PD-L1", "PD-1/PD-L1", "PD-1/PD-L1", "CTLA4", "IO+combo", "PD-1/PD-L1", "CTLA4"]')

geneSig.score <- read.csv("signature_information.csv")

# Apply a function over the loaded datasets to perform survival or response analysis
assoc.res <- lapply(1:length(expr), function(k){
    dataset <- expr[[k]]
    dataset_name <- names(expr)[k]

    if ('OS' == 'OS') {
        res.all <- lapply(1:nrow(geneSig.score), function(i) {
            sig_name <- rownames(geneSig.score)[i]
            geneSig_vector <- as.numeric(geneSig.score[i, ])
            geneSig_vector <- geneSig_vector[!is.na(geneSig_vector)]

            # Debugging output
            print(paste("Processing dataset:", dataset_name, "signature:", sig_name))

            # Check for presence of required variables
            if (!all(c("status", "time") %in% names(dataset))) {
                print(paste("Missing status or time in dataset:", dataset_name))
                return(NULL)
            }

            # Debugging output: lengths of variables
            print(paste("Lengths - status:", length(dataset$status), "time:", length(dataset$time), "variable:", length(geneSig_vector)))

            # Ensure variables have data
            if (length(dataset$status) == 0 || length(dataset$time) == 0 || length(geneSig_vector) == 0) {
                print(paste("Empty variable in dataset:", dataset_name))
                return(NULL)
            }

            # Further debug: print a few values of each variable
            print("First few values of status:")
            print(head(dataset$status))
            print("First few values of time:")
            print(head(dataset$time))
            print("First few values of geneSig_vector:")
            print(head(geneSig_vector))

            res <- geneSigSurvCont(
                dat.icb = dataset,
                geneSig = geneSig_vector,
                time.censor = 36,
                n.cutoff = 15,
                study = dataset_name,
                surv.outcome = "OS",
                sig.name = sig_name,
                cancer.type = cancer_types[k],
                treatment = treatment_types[k]
            )

            return(res)
        })
        return(res.all)
    } else if ('OS' == 'PFS') {
        res.all <- lapply(1:nrow(geneSig.score), function(i) {
            sig_name <- rownames(geneSig.score)[i]
            geneSig_vector <- as.numeric(geneSig.score[i, ])
            geneSig_vector <- geneSig_vector[!is.na(geneSig_vector)]

            # Debugging output
            print(paste("Processing dataset:", dataset_name, "signature:", sig_name))

            # Check for presence of required variables
            if (!all(c("status", "time") %in% names(dataset))) {
                print(paste("Missing status or time in dataset:", dataset_name))
                return(NULL)
            }

            # Debugging output: lengths of variables
            print(paste("Lengths - status:", length(dataset$status), "time:", length(dataset$time), "variable:", length(geneSig_vector)))

            # Ensure variables have data
            if (length(dataset$status) == 0 || length(dataset$time) == 0 || length(geneSig_vector) == 0) {
                print(paste("Empty variable in dataset:", dataset_name))
                return(NULL)
            }

            res <- geneSigSurvCont(
                dat.icb = dataset,
                geneSig = geneSig_vector,
                time.censor = 24,
                n.cutoff = 15,
                study = dataset_name,
                surv.outcome = "PFS",
                sig.name = sig_name,
                cancer.type = cancer_types[k],
                treatment = treatment_types[k]
            )

            return(res)
        })
        return(res.all)
    } else if ('OS' == 'Response') {
        res.all <- lapply(1:nrow(geneSig.score), function(i) {
            sig_name <- rownames(geneSig.score)[i]
            geneSig_vector <- as.numeric(geneSig.score[i, ])
            geneSig_vector <- geneSig_vector[!is.na(geneSig_vector)]

            # Debugging output
            print(paste("Processing dataset:", dataset_name, "signature:", sig_name))

            # Check for presence of required variables
            if (!all(c("status", "time") %in% names(dataset))) {
                print(paste("Missing status or time in dataset:", dataset_name))
                return(NULL)
            }

            # Debugging output: lengths of variables
            print(paste("Lengths - status:", length(dataset$status), "time:", length(dataset$time), "variable:", length(geneSig_vector)))

            # Ensure variables have data
            if (length(dataset$status) == 0 || length(dataset$time) == 0 || length(geneSig_vector) == 0) {
                print(paste("Empty variable in dataset:", dataset_name))
                return(NULL)
            }

            res <- geneSigLogReg(
                dat.icb = dataset,
                geneSig = geneSig_vector,
                n.cutoff = 10,
                study = dataset_name,
                sig.name = sig_name,
                n0.cutoff = 3, 
                n1.cutoff = 3,
                cancer.type = cancer_types[k],
                treatment = treatment_types[k]
            )

            return(res)
        })
        return(res.all)
    } else {
        print(paste("You can only use 'OS', 'PFS', or 'Response' as 'OS' input."))
        return(NULL)
    }
})

assoc.res <- do.call(rbind, assoc.res)
assoc.res$FDR <- p.adjust(assoc.res$Pval, method="BH")
assoc.res <- assoc.res[order(assoc.res$FDR), ]

# Meta-analysis for a gene across datasets
res_meta_pancancer <- metafun(coef = assoc.res$Coef, se = assoc.res$SE, study = assoc.res$Study, pval = assoc.res$Pval, n = assoc.res$N, cancer.type = assoc.res$Cancer_type, treatment = assoc.res$Treatment, feature = "CXCL9", cancer.spec = FALSE, treatment.spec = FALSE)

# Save the results to a CSV file
write.csv(data.frame(res_meta_pancancer), file = "Meta_Sig_analysis_OS_CXCL9_pancancer.csv", row.names = FALSE)
