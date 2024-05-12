# Use Rocker/R-ver as the base image
FROM rocker/r-ver:4.3.2

# Install system dependencies for Nextflow and PostgreSQL
RUN apt-get update && apt-get install -y \
    curl \
    unzip \
    openjdk-11-jre-headless \
    git \
    libpq-dev \
    sudo && \
    rm -rf /var/lib/apt/lists/*  # This ensures the RUN command remains a single layer

# Install Nextflow
RUN curl -s https://get.nextflow.io | bash && \
    mv nextflow /usr/local/bin/ && \
    chmod +x /usr/local/bin/nextflow


# Install R packages and skip errors
RUN install2.r --error FALSE --deps TRUE \
    readr \
    dplyr \
    meta \
    metafor \
    forestplot \
    ggplot2 \
    ggrepel \
    gridExtra \
    data.table \
    kableExtra \
    summarytools || true

# Install BiocManager
RUN R -e "install.packages('BiocManager', repos = 'http://cran.rstudio.com/')"

# Install Bioconductor packages
RUN R -e "BiocManager::install(c('GSVA', 'MultiAssayExperiment', 'survcomp', 'SummarizedExperiment'))"

# Set up environment variables
ENV PATH="/usr/local/bin:$PATH"

# Create and set permissions for R library directory
RUN mkdir -p /usr/local/lib/R/site-library && \
    chown -R root:staff /usr/local/lib/R/site-library

# Set up working directories
WORKDIR /workspace
