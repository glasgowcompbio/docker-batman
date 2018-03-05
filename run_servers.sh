#!/bin/bash

/usr/lib/rstudio-server/bin/rserver --server-app-armor-enabled=0
jupyter notebook --ip=0.0.0.0 --port=9999 --no-browser --allow-root