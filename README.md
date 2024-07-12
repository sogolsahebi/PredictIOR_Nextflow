
# PredictioR Nextflow Pipeline

## Overview
The PredictioR Nextflow pipeline is designed to analyze immunotherapy responses and identify biomarkers across various cancers. It utilizes Nextflow for workflow management and Docker for reproducibility, focusing on handling SummarizedExperiment objects for in-depth biomarker analysis.

## Software Requirements and Installation Instructions

### Nextflow
- **Version:** 24.04.2
- **Installation and Resources:**
  - **Setup Instructions:** For detailed installation steps, please refer to the [Nextflow Setup Guide](https://www.nextflow.io/docs/latest/install.html).
  - **Documentation and Resources:**
    - [General Documentation](https://www.nextflow.io/docs/latest/index.html)
    - [Nextflow Training](https://training.nextflow.io)

### Docker
- **Purpose:** Ensures computational reproducibility by containerizing the environment.
- **Installation Guide:** [Install Docker](https://docs.docker.com/get-docker/)
- **PredictioR Docker Image:**
  ```bash
  docker pull bhklab/nextflow-env
  ```
  - Docker Hub: [bhklab/nextflow-env](https://hub.docker.com/r/bhklab/nextflow-env)

## Reference Resources
- **GitHub Repository:** [PredictioR GitHub Repository](https://github.com/bhklab/PredictioR)
- **Key Publication:** [PubMed Reference, PMID: 36055464](https://pubmed.ncbi.nlm.nih.gov/36055464/)

## Data Directory Configuration

### Gene Level Analysis
- **Input Data Directory:**
  ```bash
  params.gene_data_dir = './ICB_data'
  ```
- **Example Data Files:** Includes files such as `ICB_small_Hugo.rda`, `ICB_small_Mariathasan.rda`, which are [SummarizedExperiment objects](https://bioconductor.org/packages/devel/bioc/vignettes/SummarizedExperiment/inst/doc/SummarizedExperiment.html). These files are located within the `ICB_data` directory at the [bhklab PredictioR data repository](https://github.com/bhklab/PredictioR/tree/main/data).
- **Output Data Directory:**
  ```bash
  params.out_dir = './output/gene_level_output'
  ```

### Signature Level Analysis
- **Input Data Directory:**
  ```bash
  params.signature_data_dir = './SIG_data'
  ```
- **Example Data Files:** Files like `CYT_Rooney.rda`, `EMT_Thompson.rda`, `PredictIO_Bareche.rda` are data frames with columns like:
  - `signature_name`: Name of the signature
  - `gene_name`: Name of the gene
  - `weight`: Weight assigned to each gene within the signature
  
  To see other columns, these files are also sourced from the [bhklab SignatureSets GitHub repository](https://github.com/bhklab/SignatureSets).
  The `.rda` files are stored in the object `sig` as data frames. Please follow the same format for consistency.
  
- **Output Data Directory:**
  ```bash
  params.out_dir = './output/signature_level_output'
  ```

### Meta Analysis
- **Input Data Directory:** 
  - This step aggregates the results from both gene-level and signature-level analyses.
  - **Input Directories:** 
    - Gene level output: `./output/gene_level_output`
    - Signature level output: `./output/signature_level_output`
- **Output Data Directory:**
  ```bash
  params.out_dir = './output/meta_analysis_output'
  ```

## Input Data Specifications
Ensure that clinical data is properly organized with all required and additional fields to ensure the integrity of the analysis.
- **Required Columns:**
  - `patientid`: Unique identifier for patients
  - `treatmentid`: Details of the treatment regimen
  - `response`: Patient response to treatment (Responder 'R', Non-responder 'NR')
  - `tissueid`: Standardized cancer type
  - `survival_time_pfs`: Time to progression-free survival, Example: 2.6 months
  - `survival_time_os`: Time to overall survival
  - `survival_unit`: Measurement units for survival times, typically months
  - `event_occurred_pfs`: Binary indicator of event occurrence during PFS (1,0)
  - `event_occurred_os`: Binary indicator of event occurrence during OS (1,0)

- **Additional Recommended Fields:**
  Include sex, age, histo (histological type), stage of cancer, dna, and rna details among others as necessary.

Here is the revised section focused on the signature information CSV and the related columns:

## Signature Information
For detailed information on the signatures, refer to the signature information CSV available at: [Signature Information CSV](https://github.com/bhklab/SignatureSets/tree/main/data-raw). Key columns in the CSV include:

- `signature`  Name of the signature , same names located in './SIG_data'
- `method` used for signature score calculation
- `score function`, specifying the function that should be used in the R script

For other columns and additional information, you can refer to the [Signature Information CSV](https://github.com/bhklab/SignatureSets/tree/main/data-raw).

Run the pipeline with the configured parameters using Nextflow:
```bash
nextflow run gene_level_analysis.nf
nextflow run signature_level_analysis.nf
nextflow run meta_analysis.nf
```

## Additional Notes
- Ensure that all the necessary R packages and dependencies are installed as specified in the `load_libraries.R` script.
- Customize the `nextflow.config` file to specify any additional parameters or configurations required for your specific analysis needs.







