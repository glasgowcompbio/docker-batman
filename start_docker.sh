#!/bin/bash

CODE_DIR='/Users/joewandy/git/pyBatman'
NMR_DIR='/Users/joewandy/Dropbox/Meta_clustering/NMR'

# docker run -v ${CODE_DIR}:/home/rstudio/codes -v ${NMR_DIR}:/home/rstudio/NMR -it -p 8787:8787 -p 9999:9999 batman
docker run -v ${CODE_DIR}:/home/rstudio/codes -v ${NMR_DIR}:/home/rstudio/NMR -it -p 8787:8787 -p 9999:9999 batman bash