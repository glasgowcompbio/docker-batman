docker-batman
=============

This docker image provides an automated pipeline to infer the concentration of metabolites from 1D-NMR data. It's built upon a set of Python script for the pipeline, relying on BATMAN [1] for inference.

[1] Hao, Jie, et al. "BATMAN—an R package for the automated quantification of metabolites from nuclear magnetic resonance spectra using a Bayesian model." Bioinformatics 28.15 (2012): 2088-2090

To build this image:

    $ docker build -t joewandy/docker-batman .

To install and run the image, use the following command from the shell:

    $ docker run \
    -v /Users/joewandy/git/pyBatman:/home/rstudio/NMR/codes \
    -v /Users/joewandy/Dropbox/Meta_clustering/NMR/data/test_spectra/:/home/rstudio/NMR/spectra \
    -v /Users/joewandy/Dropbox/Meta_clustering/NMR/data/background/:/home/rstudio/NMR/background \
    -v /Users/joewandy/Dropbox/Meta_clustering/NMR/results/:/home/rstudio/NMR/results \
    --name docker-batman -d -p 8787:8787 -p 9999:9999 joewandy/docker-batman

Explanation of the command above:
- `-v xxxx:yyyy` maps the host folder xxxx to yyyy in the container. You need to map the 'spectra', 'background', 'results' and 'codes' folders from your host to the container. These are where the pipeline will look for the Bruker spectra to process, the background signal to use for background correction, where to place the results and find the Python codes required to run the pipeline.
- `--name batman` gives the running container the name 'batman'.
- `-d` runs the container in a detached mode.
- `-dp 8787:8787` maps port 8787 in the host to port 8787 in the container, similarly `-dp 9999:9999` maps port 9999.
- `joewandy/docker-batman` specifies the name of the docker image to pull from dockerhub.

Please put all your spectra (in Bruker format) that you want to infer the concentrations from in the 'spectra' folder of the host. The background signals should go into the 'background' of the host. Codes can be checked out from the github repo https://github.com/joewandy/pyBatman in the host.

To attach to the shell in the running container

    $ docker exec -it docker-batman /bin/bash

And finally to access R-studio in the running container, go to `http://localhost:8787/` in the browser and log in using the username 'rstudio' and password 'rstudio'. To access Jupyter notebook in the container, go to `http://localhost:9999`.

[![](https://images.microbadger.com/badges/image/joewandy/docker-batman.svg)](https://microbadger.com/images/joewandy/docker-batman "Get your own image badge on microbadger.com")