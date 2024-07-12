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
  - This step aggregates the directories from both gene-level and signature-level analyses.
  - **Input Directories:** 
    - Gene level : `./ICB_data`
    - Signature level : `./output/signature_level_output`
- **Output Data Directory:**
  ```bash
  params.out_dir = './output/meta_analysis_output'
  ```

### Input Data Specifications
- #### ICB Data Information

  This table summarizes each dataset by study and treatment type, along with cancer types, clinical and molecular data availability, and relevant PMID 
  references. Required columns include 'treatment' and 'cancer type'.
  
  | Dataset                | Patients [#] | Cancer type | Treatment                   | Clinical endpoints | Molecular data | PMID      |
  |------------------------|--------------|-------------|-----------------------------|--------------------|----------------|-----------|
  | ICB_small_Hugo         | 27           | Melanoma    | PD-1/PD-L1                  | OS                 | RNA            | 26997480  |
  | ICB_small_Liu          | 121          | Melanoma    | PD-1/PD-L1                  | PFS/OS             | RNA/DNA        | 31792460  |
  | ICB_small_Miao         | 33           | Kidney      | PD-1/PD-L1                  | PFS/OS             | RNA/DNA        | 29301960  |
  | ICB_small_Nathanson    | 24           | Melanoma    | CTLA4                       | OS                 | RNA/DNA        | 27956380  |
  | ICB_small_Padron       | 45           | Pancreas    | PD-1/PD-L1                  | PFS/OS             | RNA            | 35662283  |
  | ICB_small_Riaz         | 46           | Melanoma    | PD-1/PD-L1                  | OS                 | RNA/DNA        | 29033130  |
  | ICB_small_Van_Allen    | 42           | Melanoma    | CTLA4                       | PFS/OS             | RNA/DNA        | 26359337  |
  | ICB_small_Mariathasan  | 195          | Bladder     | PD-1/PD-L1                  | OS                 | RNA/DNA        | 29443960  |
  
  
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
  

- #### Signature Information

  This table summarizes each signature name by study and PMID references, the method for computing the signature score, and the corresponding score 
  function.
  
  | Signature            | DNA/RNA | RNA Type           | Method | Cancer Type         | Score Function | PMID      |
  |----------------------|---------|--------------------|--------|---------------------|----------------|-----------|
  | ADO_Sidders          | RNA     | Count RNA-seq/TPM  | GSVA   | Multiple            | geneSigGSVA    | 31953314  |
  | APM_Thompson         | RNA     | log CPM            | GSVA   | Lung, melanoma      | geneSigGSVA    | 33028693  |
  | APM_Wang             | RNA     | Microarray         | GSVA   | Multiple            | geneSigGSVA    | 31767055  |
  | Bcell_Budczies       | RNA     | Microarray         | GSVA   | Lung                | geneSigGSVA    | 33520406  |
  | Bcell_Helmink        | RNA     | log FPKM           | GSVA   | Melanoma, kidney    | geneSigGSVA    | 31942075  |
  | Blood_Friedlander    | RNA     | Microarray         | GSVA   | Melanoma            | geneSigGSVA    | 28807052  |
  | C-ECM_Chakravarthy   | RNA     | Normalized counts  | ssGSEA | Multiple            | geneSigssGSEA  | 30410077  |
  | CCL5-CXCL9_Dangaj    | RNA     |                    | GSVA   | Multiple            | geneSigGSVA    | 31185212  |
  | CD39-CD8Tcell_Chow   | RNA     | RNA-seq count      | GSVA   | Lung                | geneSigGSVA    | 36574773  |


  **Required Columns:**
  - `signature`: Name of the signature, same names located in './SIG_data'
  - `method`: Used for signature score calculation
  - `score function`: Specifying the function that should be used in the R script
  For detailed information on the signatures used in the pipeline, refer to the signature(there are more than 50) information CSV available at: [Signature Information CSV](https://github.com/bhklab/SignatureSets/tree/main/data-raw). 

## Running the Pipeline
Run the pipeline with the configured parameters using Nextflow:
```bash
nextflow run gene_level_analysis.nf
nextflow run signature_level_analysis.nf
nextflow run meta_analysis.nf
```

## Additional Notes
- Necessary R packages and dependencies are installed as specified in `load_libraries.R` and included in the BHK Docker.
- Customize the nextflow.config file to specify any additional parameters or configurations required for your specific analysis needs
