#!/bin/bash

CODE_DIR="/home/rstudio/codes/notebooks"
OUTPUT_FILE="/home/rstudio/NMR/results/diagnostic.html"
TEMP_DIR="/home/rstudio/NMR/temp"

/usr/games/cowsay 'Starting analysis'

jupyter nbconvert --to HTML --ExecutePreprocessor.timeout=-1 --execute ${CODE_DIR}/predict_spectra_concentrations.ipynb --output ${OUTPUT_FILE} && rm -rf ${TEMP_DIR} && figlet 'FINISH'