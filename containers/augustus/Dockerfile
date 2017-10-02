FROM ubuntu:latest

MAINTAINER Anton Kulaga <antonkulaga@gmail.com>

RUN apt-get update && apt-get install -qqy build-essential

ADD http://bioinf.uni-greifswald.de/augustus/binaries/augustus.current.tar.gz /opt/
RUN cd /opt && \
    tar -xzvf augustus* && \
    rm -rf *.tar.gz && \
    cd augustus && \
    make

ADD fgram_base /opt/augustus/config/species/
ENV AUGUSTUS_CONFIG_PATH /opt/augustus/config
ENV PATH /usr/local/sbin:/usr/local/bin:/usr/sbin:/usr/bin:/sbin:/bin:/opt/augustus/bin:/opt/augustus/scripts
CMD["/opt/augustus/bin/augustus"]