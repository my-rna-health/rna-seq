#!/usr/bin/env bash
docker run -v `pwd`/data:/home/rstudio/data/ -p 8787:8787 sleuth:latest
