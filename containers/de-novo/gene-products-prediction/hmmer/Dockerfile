FROM quay.io/comp-bio-aging/augustus:latest

MAINTAINER Anton Kulaga

ENV HMMER_VERSION="3.1b2"
ENV HMMER_PLATFORM="linux-intel-x86_64"

WORKDIR /opt/

RUN wget http://eddylab.org/software/hmmer3/3.1b2/hmmer-$HMMER_VERSION-$HMMER_PLATFORM.tar.gz && \
    tar -zxvf hmmer-$HMMER_VERSION-$HMMER_PLATFORM.tar.gz && \
    rm -rf hmmer-$HMMER_VERSION-$HMMER_PLATFORM.tar.gz && \
    ln -s hmmer-$HMMER_VERSION-$HMMER_PLATFORM hmmer

WORKDIR /opt/hmmer

ENV PATH "$PATH:/opt/hmmer:/opt/hmmer/binaries"

RUN mkdir /data
WORKDIR /data

CMD ["/bin/bash"]