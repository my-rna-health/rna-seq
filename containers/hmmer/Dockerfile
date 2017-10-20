FROM augustus:latest

MAINTAINER Anton Kulaga

ENV HMMER_VERSION="3.1b2"
ENV HMMER_PLATFORM="linux-intel-x86_64"

RUN apt-get update && \
 apt-get upgrade -y && \
 apt-get install -y cmake zip wget gcc git build-essential

WORKDIR /opt/

RUN wget http://eddylab.org/software/hmmer3/3.1b2/hmmer-$HMMER_VERSION-$HMMER_PLATFORM.tar.gz && \
    tar -zxvf hmmer-$HMMER_VERSION-$HMMER_PLATFORM.tar.gz && \
    rm -rf hmmer-$HMMER_VERSION-$HMMER_PLATFORM.tar.gz && \
    ln -s hmmer-$HMMER_VERSION-$HMMER_PLATFORM hmmer

WORKDIR /opt/hmmer

ENV PATH "$PATH:/opt/hmmer"

RUN mkdir /data
WORKDIR /data

CMD ["/bin/bash"]