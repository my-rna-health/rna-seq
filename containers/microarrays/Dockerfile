FROM rocker/rstudio

MAINTAINER antonkulaga@gmail.com

# nuke cache dirs before installing pkgs; tip from Dirk E fixes broken img
RUN  rm -f /var/lib/dpkg/available && rm -rf  /var/cache/apt/*

RUN apt-get update && \
    apt-get -y  install --fix-missing gdb libxml2-dev python-pip libz-dev software-properties-common python3-software-properties automake nano cmake zip wget gcc git build-essential curl libmariadb-client-lgpl-dev

RUN pip install awscli

ADD install.R /tmp/

RUN R -f /tmp/install.R && \
    echo "library(BiocInstaller)" > $HOME/.Rprofile

# invalidates cache every 24 hours
ADD http://master.bioconductor.org/todays-date /tmp/