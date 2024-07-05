
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
  docker pull sogolsahebi/nextflow-env
  ```
  - Docker Hub: [sogolsahebi/nextflow-env](https://hub.docker.com/r/sogolsahebi/nextflow-env)

## Reference Resources
- **GitHub Repository:** [PredictioR GitHub Repository](https://github.com/bhklab/PredictioR)
- **Key Publication:** [PubMed Reference, PMID: 36055464](https://pubmed.ncbi.nlm.nih.gov/36055464/)

## Data Directory Configuration

### Signature Level Analysis
- **Input Data Directory:**
  ```bash
  params.signature_data_dir = './ICB_data'
  ```
- **Example Data Files:** Includes files such as `ICB_small_Hugo.rda`, `ICB_small_Mariathasan.rda`, which are [SummarizedExperiment objects](https://bioconductor.org/packages/devel/bioc/vignettes/SummarizedExperiment/inst/doc/SummarizedExperiment.html). These files are located within the `ICB_data` directory at the [bhklab PredictioR data repository](https://github.com/bhklab/PredictioR/tree/main/data).
- **Output Data Directory:**
  ```bash
  params.signature_out_dir = './output'
  ```

### Gene Level Analysis
- **Input Data Directory:**
  ```bash
  params.gene_data_dir = './SIG_data'
  ```
- **Example Data Files:** Files like `CYT_Rooney.rda`, `EMT_Thompson.rda`, `PredictIO_Bareche.rda` are data frames with columns:
  - `signature_name`: Name of the signature
  - `gene_id`: Gene identifier
  - `gene_type`: Type of gene
  - `gene_name`: Name of the gene
  - `weight`: Weight assigned to each gene within the signature
    
  These files are also sourced from the [bhklab SignatureSets GitHub repository](https://github.com/bhklab/SignatureSets).
- **Output Data Directory:**
  ```bash
  params.gene_out_dir = './output'
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

## Execution Instructions
Run the pipeline with the configured parameters using Nextflow:
```bash
nextflow run signature_level_analysis.nf
nextflow run gene_level_analysis.nf
```

