# Dockerfile to install rstudio-server and BATMAN library for NMR spectral processing
FROM ubuntu:14.04
MAINTAINER Joe Wandy <joe.wandy@glasgow.ac.uk>

RUN apt-get -y update \
      && apt-get -y upgrade \
      && apt-get -y install curl wget r-base libapparmor1 libcurl4-openssl-dev libxml2-dev libssl-dev gdebi-core \
      && apt-cache search r-cran | cut -f 1 -d ' ' | xargs apt-get -y install

# grab latest rstudio-server
RUN curl https://s3.amazonaws.com/rstudio-server/current.ver | \
        xargs -I {} wget http://download2.rstudio.org/rstudio-server-{}-amd64.deb -O rstudio.deb \
      && gdebi -n rstudio.deb \
      && rm rstudio.deb \
      && apt-get clean

# install ssh, required for Rmpi (parallelisation)
RUN apt-get -y install ssh

# better top
RUN apt-get -y install htop
ENV TERM xterm

# install batman & other R packages here
ADD installBatman.R /home/root/installBatman.R
RUN /usr/bin/Rscript /home/root/installBatman.R

# create rstudio user and copy the example folder into the container
RUN useradd -m -d /home/rstudio rstudio \
      && echo rstudio:rstudio | chpasswd
ADD example /home/rstudio/example
RUN chown -R rstudio:rstudio /home/rstudio/example

# expose port 8787 to the host for rstudio-server
EXPOSE 8787

# start rstudio-server in the background
CMD ["/usr/lib/rstudio-server/bin/rserver", "--server-daemonize=0", "--server-app-armor-enabled=0"]
