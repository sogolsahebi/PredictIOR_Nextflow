# Use Rocker/R-ver as the base image
FROM rocker/r-ver:4.3.2

<<<<<<< HEAD
=======

>>>>>>> b893ab0 (Repo is ready to use)
# Install system dependencies for Nextflow, Docker, and PostgreSQL
RUN apt-get update && apt-get install -y \
    curl \
    unzip \
    openjdk-11-jre-headless \
    git \
    libpq-dev \
    sudo \
    bash \
    apt-transport-https \
    ca-certificates \
    software-properties-common && \
    rm -rf /var/lib/apt/lists/*

# Install Docker
RUN curl -fsSL https://download.docker.com/linux/ubuntu/gpg | sudo apt-key add - && \
    sudo add-apt-repository "deb [arch=amd64] https://download.docker.com/linux/ubuntu $(lsb_release -cs) stable" && \
    sudo apt-get update && \
    sudo apt-get install -y docker-ce && \
    sudo usermod -aG docker $(whoami) && \
    rm -rf /var/lib/apt/lists/*

# Install Nextflow
RUN curl -s https://get.nextflow.io | bash && \
    mv nextflow /usr/local/bin/ && \
    chmod +x /usr/local/bin/nextflow

# Verify Nextflow installation
RUN ls -l /usr/local/bin/nextflow && which nextflow

# Install necessary system dependencies for R packages
RUN apt-get update && apt-get install -y \
    libcurl4-gnutls-dev \
    libssl-dev \
    libxml2-dev \
    libblas-dev \
    liblapack-dev \
    libbz2-dev \
    liblzma-dev \
    zlib1g-dev && \
    rm -rf /var/lib/apt/lists/*

# Install Bioconductor and CRAN packages
RUN R -e "install.packages('BiocManager', repos='http://cran.rstudio.com/')"
RUN R -e "BiocManager::install(c('SummarizedExperiment', 'MultiAssayExperiment', 'GSVA', 'survcomp'))"
RUN R -e "install.packages(c('readr', 'dplyr', 'meta', 'metafor', 'forestplot', 'ggplot2', 'ggrepel', 'gridExtra', 'data.table', 'kableExtra', 'survival', 'prodlim'), dependencies=TRUE)"
RUN R -e "install.packages('box', repos = 'https://klmr.r-universe.dev')"

# Verify R installation
RUN which R

# Set up environment variables
ENV PATH="/usr/local/bin:$PATH"

# Create and set permissions for R library directory
RUN mkdir -p /usr/local/lib/R/site-library && \
    chown -R root:staff /usr/local/lib/R/site-library

# Set up the working directory
WORKDIR /PredictIOR_Nextflow
