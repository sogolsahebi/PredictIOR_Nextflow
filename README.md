# PredictioR Nextflow Pipeline

## Overview
The PredictioR Nextflow pipeline is specifically designed to analyze immunotherapy responses and identify biomarkers across various cancers. This comprehensive guide covers all aspects of the pipeline's setup and operation.

## Software Requirements and Installation Instructions

### Nextflow
- **Version:** 24.04.2
- **Installation:**
  Follow the official [Nextflow Installation Guide](https://www.nextflow.io/docs/latest/install.html) for detailed steps.
- **Documentation and Resources:**
  - [General Documentation](https://www.nextflow.io/docs/latest/index.html)
  - [Mentorship Programs](https://nf-co.re/mentorships)
  - [Training Opportunities](https://training.nextflow.io)

### Docker
- **Purpose:** Ensures computational reproducibility by containerizing the environment.
- **Installation Guide:** [Install Docker](https://docs.docker.com/get-docker/)
- **PredictioR Docker Image:**
  ```bash
  docker pull sogolsahebi/nextflow-env
  ```
  Docker Hub: [sogolsahebi/nextflow-env](https://hub.docker.com/r/sogolsahebi/nextflow-env)

## Reference Resources
- **GitHub Repository:** [PredictioR GitHub Repository](https://github.com/bhklab/PredictioR)
- **Key Publication:** [PubMed Reference PMID: 36055464](https://pubmed.ncbi.nlm.nih.gov/36055464/)

## Data Directory Configuration

### Signature Level Analysis
- **Input Data Directory:**
  ```bash
  params.signature_data_dir = './ICB_data'
  ```
  Example Data Files: Includes `ICB_small_Hugo.rda`, `ICB_small_Mariathasan.rda`, etc.
- **Output Data Directory:**
  ```bash
  params.signature_out_dir = './output'
  ```

### Gene Level Analysis
- **Input Data Directory:**
  ```bash
  params.gene_data_dir = './SIG_data'
  ```
  Example Data Files: Includes `CYT_Rooney.rda`, `EMT_Thompson.rda`, `PredictIO_Bareche.rda`, etc.
- **Output Data Directory:**
  ```bash
  params.gene_out_dir = './output'
  ```

## Configuration of Pipeline Parameters

Configure the pipeline based on the specifics of the cancer type and treatment modalities:

- **Cancer Type Configuration:**
  ```bash
  cancer_type = 'Bladder'  // Specify the cancer type for analysis
  ```

- **Treatment Type Configuration:**
  ```bash
  treatment = 'PD-1/PD-L1'  // Define the treatment regimen
  ```

## Input Data Specifications

Ensure clinical data is organized accurately to maintain the integrity of the analysis:

- **Required Columns:**
  - `patientid`: Unique identifier for patients.
  - `treatmentid`: Details of the treatment regimen.
  - `response`: Patient response to treatment (Responder 'R', Non-responder 'NR').
  - `tissueid`: Standardized cancer type.
  - `survival_time_pfs`: Time to progression-free survival.
  - `survival_time_os`: Time to overall survival.
  - `survival_unit`: Measurement units for survival times, typically months.
  - `event_occurred_pfs`: Binary indicator of event occurrence during PFS.
  - `event_occurred_os`: Binary indicator of event occurrence during OS.

- **Additional Recommended Fields:**
  Include `sex`, `age`, `histo` (histological type), `stage` of cancer, `dna`, and `rna` details among others as necessary.

## Execution Instructions

Run the pipeline with the configured parameters using Nextflow:

```bash
nextflow run signature_level_analysis.nf
nextflow run gene_level_analysis.nf
```

