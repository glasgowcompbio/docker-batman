#!/bin/bash

/usr/lib/rstudio-server/bin/rserver --server-app-armor-enabled=0
jupyter notebook --port=9999 --no-browser --ip=0.0.0.0