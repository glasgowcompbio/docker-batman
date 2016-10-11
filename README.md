docker-batman
=============

This docker container provides BATMAN [1] running on top of Ubuntu 14.04 and R-Studio, used for the Bayesian analysis of 1D-NMR data. It seems a bit difficult to get BATMAN to compile on the latest OSX version, so I created this container ..

[1] Hao, Jie, et al. "BATMAN—an R package for the automated quantification of metabolites from nuclear magnetic resonance spectra using a Bayesian model." Bioinformatics 28.15 (2012): 2088-2090

To install and run the docker BATMAN container:

    $ docker run -v /Users/joewandy/Dropbox/Analysis/NMR:/home/rstudio/NMR --name batman -dp 8787:8787 joewandy/docker-batman

Explanations:
- This is to map the host folder /Users/joewandy/Dropbox/Analysis/NMR to /home/rstudio/NMR in the container
-v /Users/joewandy/Dropbox/Analysis/NMR:/home/rstudio/NMR
- This is to give the running container a name 'batman'
--name batman
- Maps port 8787 in the host to port 8787 on the container
-dp 8787:8787
- Specify the name of the docker image hosted in the dockerhub to pull
joewandy/docker-batman

To go into the shell of the running batman container

    $ docker exec -it batman /bin/bash

[![](https://images.microbadger.com/badges/image/joewandy/docker-batman.svg)](https://microbadger.com/images/joewandy/docker-batman "Get your own image badge on microbadger.com")