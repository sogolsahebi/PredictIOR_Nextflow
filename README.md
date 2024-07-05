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
  - [Docker Hub: sogolsahebi/nextflow-env](https://hub.docker.com/r/sogolsahebi/nextflow-env)

## Reference Resources
- **GitHub Repository:** [PredictioR GitHub Repository](https://github.com/bhklab/PredictioR)
- **Key Publication:** [PubMed Reference, PMID: 36055464](https://pubmed.ncbi.nlm.nih.gov/36055464/)

## Data Directory Configuration

### Gene Level Analysis
- **Input Data Directory:**
  ```bash
  params.gene_data_dir = './SIG_data'
  ```
- **Example Data Files:** Includes `ICB_small_Hugo.rda`, `ICB_small_Mariathasan.rda`, which are [**SummarizedExperiment** objects](https://bioconductor.org/packages/devel/bioc/vignettes/SummarizedExperiment/inst/doc/SummarizedExperiment.html#anatomy-of-a-summarizedexperiment). These files are located in the [bhklab PredictioR data repository](https://github.com/bhklab/PredictioR/tree/main/data).

- **Output Data Directory:**
  ```bash
  params.gene_out_dir = './output'
  ```

### Signature Level Analysis
- **Input Data Directory:**
  ```bash
  params.signature_data_dir = './ICB_data'
  ```
- **Example Data Files:** Files such as `CYT_Rooney.rda`, `EMT_Thompson.rda`, and `PredictIO_Bareche.rda` are located in the [bhklab SignatureSets GitHub repository](https://github.com/bhklab/SignatureSets). These files are data frames containing information about our signatures, structured with the following columns:
  - `signature_name`: The name of the signature.
  - `gene_id`: Identifier for the gene.
  - `gene_type`: Type of the gene.
  - `gene_name`: Name of the gene.
  - `weight`: Weight assigned to the gene within the signature.

- **Output Data Directory:**
  ```bash
  params.signature_out_dir = './output'
  ```

## Handling SummarizedExperiment Objects

SummarizedExperiment objects are integral for our analyses, particularly in the signature level pipeline:
- **Starting Example:** `ICB_small_Mariathasan.rda` used in `signature_level_analysis.nf`.
- [Bioconductor SummarizedExperiment Documentation](https://bioconductor.org/packages/devel/bioc/vignettes/SummarizedExperiment/inst/doc/SummarizedExperiment.html).

## Configuration of Pipeline Parameters
Configure the pipeline for specific cancer types and treatments:
- **Cancer Type:**
  ```bash
  cancer_type = 'Bladder'  // Specify the cancer type for analysis
  ```
- **Treatment Type:**
  ```bash
  treatment = 'PD-1/PD-L1'  // Define the treatment regimen
  ```

## Input Data Specifications
Organize clinical data with required and additional fields to ensure the integrity of the analysis.

## Execution Instructions
Run the pipeline with the configured parameters using Nextflow:
```bash
nextflow run signature_level_analysis.nf
nextflow run gene_level_analysis.nf

