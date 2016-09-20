docker-batman
=============

This docker container provides BATMAN [1] running on top of Ubuntu 14.04 and R-Studio, used for the Bayesian analysis of 1D-NMR data. It seems a bit difficult to get BATMAN to compile on the latest OSX version, so I created this container ..

[1] Hao, Jie, et al. "BATMAN—an R package for the automated quantification of metabolites from nuclear magnetic resonance spectra using a Bayesian model." Bioinformatics 28.15 (2012): 2088-2090

To pull this container (joewandy/docker-batman) from the hub, assign the name 'batman' to the running container and map port 8787 from the container to the host.

    $ docker run --name batman -dp 8787:8787 joewandy/docker-batman

Next, to find out the IP address assigned to your virtual machine

    $ docker-machine ip default

Using that IP, run R-Studio from the browser, e.g. http://192.168.99.100:8787. Login using user rstudio, password rstudio. An example script runBatman.R is provided inside /home/rstudio/example. Open that script and run it to test BATMAN on some serum data. The results can be found in the /home/rstudio/example/runBATMAN/BatmanOutput.

To access the shell once the batman container is running

    $ docker exec -it batman /bin/bash

[![](https://images.microbadger.com/badges/image/joewandy/docker-batman.svg)](https://microbadger.com/images/joewandy/docker-batman "Get your own image badge on microbadger.com")