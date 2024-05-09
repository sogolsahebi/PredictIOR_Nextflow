# Example Dockerfile for R environment
FROM r-base:latest
RUN R -e "install.packages('ggplot2', repos='http://cran.rstudio.com/')"

# git add .gitpod.yml Dockerfile
# git commit -m "Add Gitpod and Docker configurations"
# git push
