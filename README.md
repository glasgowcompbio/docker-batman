docker-batman
=============

This docker image provides an automated pipeline to infer the concentration of metabolites from 1D-NMR data. It's built upon a set of Python script for the pipeline, relying on BATMAN [1] for inference.

[1] Hao, Jie, et al. "BATMAN—an R package for the automated quantification of metabolites from nuclear magnetic resonance spectra using a Bayesian model." Bioinformatics 28.15 (2012): 2088-2090

To build this image:

    $ docker build -t joewandy/docker-batman .

To run the image, refer to start_docker.sh. In short:

    NMR_DIR='/Users/joewandy/Dropbox/Meta_clustering/NMR/test'
    docker run --rm -v ${NMR_DIR}:/home/rstudio/NMR --name docker-batman joewandy/docker-batman:latest

The analysis script, found at https://github.com/joewandy/pyBatman, will specifically look at /home/rstudio/NMR (in the container) to find all the input files it needs. In particular, put all your spectra (in Bruker format) in the 'spectra' folder of ${NMR_DIR}. The background signals should go into the 'background' folder in ${NMR_DIR}.

To start a bash shell in the container, run the start_docker_shell.sh script. You can access R-studio in the running container by going to `http://localhost:8787/` in the browser and log in using the username 'rstudio' and password 'rstudio'. To access Jupyter notebook in the container, go to `http://localhost:9999`.

[![](https://images.microbadger.com/badges/image/joewandy/docker-batman.svg)](https://microbadger.com/images/joewandy/docker-batman "Get your own image badge on microbadger.com")