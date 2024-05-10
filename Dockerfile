# Start from a base image with necessary tools
FROM ubuntu:20.04

# Avoid prompts from apt by setting this environment variable
ENV DEBIAN_FRONTEND=noninteractive

# Install Nextflow
# RUN curl -s https://get.nextflow.io | bash
#RUN mv nextflow /usr/local/bin/

# Install system dependencies for R, R packages, Git, and sudo
RUN apt-get update && apt-get install -y \
    software-properties-common \
    dirmngr \
    gpg-agent \
    git \
    libcurl4-gnutls-dev \
    libssl-dev \
    libxml2-dev \
    libfontconfig1-dev \
    libcairo2-dev \
    libxt-dev \
    sudo \
    && apt-key adv --keyserver keyserver.ubuntu.com --recv-keys 0x51716619e084dab9 \
    && add-apt-repository "deb http://cloud.r-project.org/bin/linux/ubuntu $(lsb_release -cs)-cran40/" \
    && apt-get update && apt-get install -y r-base

# Clean up to reduce the image size
RUN apt-get clean && rm -rf /var/lib/apt/lists/*

# Allow passwordless sudo for all users
RUN echo "ALL ALL=(ALL) NOPASSWD: ALL" > /etc/sudoers.d/nopasswd

# Set environment variables
ENV R_HOME /usr/lib/R

# Install additional R packages
RUN R -e "install.packages(c('survcomp', 'GSVA', 'dplyr', 'meta', 'metafor', 'forestplot', 'ggplot2', 'ggrepel', 'gridExtra', 'data.table', 'kableExtra', 'summarytools', 'MultiAssayExperiment'), dependencies=TRUE, repos='http://cran.rstudio.com/')"
