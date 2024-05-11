
# Start from a Debian base image
FROM debian:buster

# Install system dependencies for R
RUN apt-get update && apt-get install -y \
    r-base \
    r-base-dev

# Verify R installation
RUN R --version

# Start from the official R base image
FROM r-base:4.3.2

#curl -s https://get.sdkman.io | bash
#sudo mv nextflow /usr/local/bin

# Install system dependencies for R packages, curl, git, and Java (for Nextflow)
RUN apt-get update && apt-get install -y \
    curl \
    git \
    libcurl4-openssl-dev \
    libssl-dev \
    libxml2-dev \
    libfontconfig1-dev \
    libcairo2-dev \
    libxt-dev \
    libbz2-dev \
    liblzma-dev \
    openjdk-11-jdk  # Installing Java

# Install Nextflow
RUN curl -s https://get.nextflow.io | bash && \
    mv nextflow /usr/local/bin/

# Copy the install_packages.R script from your R folder to the /tmp directory in the image
COPY R/install_packages.R /tmp/
