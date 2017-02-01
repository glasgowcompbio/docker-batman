#!/bin/bash

NMR_DIR='/Users/joewandy/Dropbox/Meta_clustering/NMR/test'

docker run --rm -v ${NMR_DIR}:/home/rstudio/NMR --name docker-batman joewandy/docker-batman:latest