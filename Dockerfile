# Use Rocker/R-ver as the base image that includes RStudio
FROM rocker/rstudio:4.3.2

# Disable authentication for RStudio
ENV DISABLE_AUTH=true

# Install system dependencies for Nextflow, Docker, PostgreSQL, and Pandoc
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
    software-properties-common \
    pandoc && \
    rm -rf /var/lib/apt/lists/*

# Install Docker
RUN curl -fsSL https://download.docker.com/linux/ubuntu/gpg | apt-key add - && \
    add-apt-repository "deb [arch=amd64] https://download.docker.com/linux/ubuntu $(lsb_release -cs) stable" && \
    apt-get update && \
    apt-get install -y docker-ce && \
    usermod -aG docker rstudio && \
    rm -rf /var/lib/apt/lists/*

# Install Nextflow
RUN mkdir -p /usr/local/bin && \
    if [ ! -f /usr/local/bin/nextflow ]; then \
        curl -s https://get.nextflow.io | bash && \
        mv nextflow /usr/local/bin/ && \
        chmod 755 /usr/local/bin/nextflow; \
    fi

# Install necessary R packages
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
RUN R -e "install.packages('BiocManager', repos='http://cran.rstudio.com/')" && \
    R -e "BiocManager::install(c('SummarizedExperiment', 'MultiAssayExperiment', 'GSVA', 'survcomp', 'BiocStyle', 'EnhancedVolcano'))" && \
    R -e "install.packages(c('readr', 'dplyr', 'meta', 'metafor', 'forestplot', 'ggplot2', 'ggrepel', 'gridExtra', 'data.table', 'kableExtra', 'survival', 'prodlim', 'summarytools'), dependencies=TRUE)" && \
    R -e "install.packages('box', repos = 'https://klmr.r-universe.dev')" && \
    R -e "install.packages('rmarkdown', repos='http://cran.rstudio.com/')"

# Install additional R packages
RUN R -e "install.packages(c('stringr', 'rstudioapi', 'pheatmap', 'RColorBrewer'), dependencies=TRUE)"

# Install missing CRAN and Bioconductor packages
RUN R -e "BiocManager::install(c('ComplexHeatmap', 'survminer', 'cowplot', 'ggsignif'))" && \
    R -e "install.packages(c('EnhancedVolcano'), dependencies=TRUE)"

# Install the PredictioR package from GitHub
RUN R -e "devtools::install_github('bhklab/PredictioR')"

# Set up the working directory
WORKDIR /PredictIOR_Nextflow

# Add a script to render Rmd files
COPY render_rmd.sh /usr/local/bin/render_rmd.sh
RUN chmod +x /usr/local/bin/render_rmd.sh

# Expose the RStudio port
EXPOSE 8787

# Command to run when the container starts
CMD ["/usr/lib/rstudio-server/bin/rserver", "--server-daemonize=0"]

# Ensure the PATH includes the directory where Nextflow is installed
RUN echo 'export PATH="/usr/local/bin:$PATH"' >> /etc/bash.bashrc

# Add library calling file
COPY load_libraries.R /usr/local/lib/R/site-library/load_libraries.R
RUN mkdir -p /R && ln -s /usr/local/lib/R/site-library/load_libraries.R /R/load_libraries.R
