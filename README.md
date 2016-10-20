docker-batman
=============

This docker image provides BATMAN [1] running on top of Ubuntu 16.04 and R-Studio, used for the Bayesian analysis of 1D-NMR data. It seems a bit difficult to get BATMAN to compile on the latest OSX version, so I created this container .. Also included in this image is Anaconda Python, alongside libraries commonly used for data analysis (numpy, etc.).

[1] Hao, Jie, et al. "BATMAN—an R package for the automated quantification of metabolites from nuclear magnetic resonance spectra using a Bayesian model." Bioinformatics 28.15 (2012): 2088-2090

To build this image:

    $ docker build -t joewandy/docker-batman .

To install and run the image, use the following command from the shell:

    $ docker run -v /Users/joewandy/Dropbox/Analysis/NMR:/home/rstudio/NMR --name batman -d -p 8787:8787 -p 8888:8888 joewandy/docker-batman

Explanation of the command above:
- `-v /Users/joewandy/Dropbox/Analysis/NMR:/home/rstudio/NMR` maps the host folder /Users/joewandy/Dropbox/Analysis/NMR to /home/rstudio/NMR in the container.
- `--name batman` gives the running container the name 'batman'.
- `-d` runs the container in a detached mode.
- `-dp 8787:8787` maps port 8787 in the host to port 8787 in the container, similarly `-dp 8787:8787` maps port 8888.
- `joewandy/docker-batman` specifies the name of the docker image to pull from dockerhub.

To attach to the shell in the running container

    $ docker exec -it batman /bin/bash

And finally to access R-studio in the running container, go to `http://localhost:8787/` in the browser and log in using the username 'rstudio' and password 'rstudio'. To access Jupyter notebook in the container, go to `http://localhost:8888`.

[![](https://images.microbadger.com/badges/image/joewandy/docker-batman.svg)](https://microbadger.com/images/joewandy/docker-batman "Get your own image badge on microbadger.com")