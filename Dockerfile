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

# Verify Nextflow installation
RUN which nextflow
# sudo mv nextflow /usr/local/bin


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

# Install necessary system dependencies
RUN apt-get update && apt-get install -y \
    libcurl4-gnutls-dev \
    libssl-dev \
    libxml2-dev \
    libblas-dev \
    liblapack-dev \
    libbz2-dev \
    liblzma-dev \
    zlib1g-dev

# Install Bioconductor packages
RUN R -e "install.packages('BiocManager', dependencies=TRUE)"
RUN R -e "BiocManager::install('SummarizedExperiment')"

# RUN R -e "install.packages('BiocManager'); BiocManager::install(c('GSVA', 'MultiAssayExperiment', 'survcomp', 'SummarizedExperiment'))"

# Install the 'box' package from R-universe
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

# Copy R scripts into the Docker image
COPY ./R /PredictIOR_Nextflow/R

# Expose any necessary ports, for example, if your app has a web server
# EXPOSE 8080

# Command to run on container start
CMD ["R"]

