FROM quay.io/comp-bio-aging/basic:latest

MAINTAINER Anton Kulaga

ENV HHSUITE_VERSION="3.0-beta.3"
ENV HHSUITE_PLATFORM="linux"

WORKDIR /opt/

RUN wget https://github.com/soedinglab/hh-suite/releases/download/v$HHSUITE_VERSION/hhsuite-$HHSUITE_VERSION-Linux.deb && \
    tar -zxvf hmmer-$HMMER_VERSION-$HMMER_PLATFORM.tar.gz && \
    rm -rf hmmer-$HMMER_VERSION-$HMMER_PLATFORM.tar.gz && \
    ln -s hmmer-$HMMER_VERSION-$HMMER_PLATFORM hmmer

WORKDIR /opt/hmmer

ENV PATH "$PATH:/opt/hmmer:/opt/hmmer/binaries"

RUN mkdir /data
WORKDIR /data

CMD ["/bin/bash"]