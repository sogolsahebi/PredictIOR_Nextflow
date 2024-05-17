#!/bin/bash

# Check if an argument is provided
if [ -z "$1" ]; then
  echo "No Rmd file provided. Usage: render_rmd.sh <path_to_Rmd_file>"
  exit 1
fi

RMD_FILE=$1

# Render the Rmd file
R -e "rmarkdown::render('$RMD_FILE')"
