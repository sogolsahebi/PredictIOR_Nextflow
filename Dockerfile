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
RUN curl -fsSL https://download.docker.com/linux/ubuntu/gpg | sudo apt-key add - && \
    sudo add-apt-repository "deb [arch=amd64] https://download.docker.com/linux/ubuntu $(lsb_release -cs) stable" && \
    sudo apt-get update && \
    sudo apt-get install -y docker-ce && \
    sudo usermod -aG docker rstudio && \
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
    R -e "BiocManager::install(c('SummarizedExperiment', 'MultiAssayExperiment', 'GSVA', 'survcomp', 'BiocStyle'))" && \
    R -e "install.packages(c('readr', 'dplyr', 'meta', 'metafor', 'forestplot', 'ggplot2', 'ggrepel', 'gridExtra', 'data.table', 'kableExtra', 'survival', 'prodlim'), dependencies=TRUE)" && \
    R -e "install.packages('box', repos = 'https://klmr.r-universe.dev')" && \
    R -e "install.packages('rmarkdown', repos='http://cran.rstudio.com/')"

# Set up the working directory
WORKDIR /PredictIOR_Nextflow

# Add a script to render Rmd files
COPY render_rmd.sh /usr/local/bin/render_rmd.sh
RUN chmod +x /usr/local/bin/render_rmd.sh

# Expose the RStudio port
EXPOSE 8787

# Command to run when the container starts
CMD ["/usr/lib/rstudio-server/bin/rserver", "--server-daemonize=0"]
