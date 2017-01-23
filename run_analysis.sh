#!/bin/bash

jupyter nbconvert --to HTML --ExecutePreprocessor.timeout=-1 --execute /home/rstudio/codes/notebooks/predict_spectra_concentrations.ipynb