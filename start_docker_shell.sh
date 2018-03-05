#!/bin/bash

CODE_DIR='/Users/joewandy/git/pyBatman'
NMR_DIR='/Users/joewandy/Dropbox/Analysis/NMR/promastigotes/'

docker run --rm -v ${CODE_DIR}:/home/rstudio/codes -v ${NMR_DIR}:/home/rstudio/NMR --name docker-batman -it -p 8787:8787 -p 9999:9999 joewandy/docker-batman:latest bash
