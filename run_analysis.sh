#!/bin/bash

CODE_DIR="/home/rstudio/codes/notebooks"

/usr/games/cowsay 'Starting analysis'
echo
python ${CODE_DIR}/start_analysis.py && figlet 'FINISH'