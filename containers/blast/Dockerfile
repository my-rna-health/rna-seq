FROM quay.io/comp-bio-aging/hmmer:latest

MAINTAINER Anton Kulaga

ENV BLAST_VERSION="2.6.0"
ENV BLAST_PLATFORM="x64-linux"

WORKDIR /opt/

RUN wget ftp://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/${BLAST_VERSION}/ncbi-blast-${BLAST_VERSION}+-${BLAST_PLATFORM}.tar.gz && \
    tar -zxvf ncbi-blast-${BLAST_VERSION}+-${BLAST_PLATFORM}.tar.gz && \
    rm -rf ncbi-blast-${BLAST_VERSION}+-${BLAST_PLATFORM}.tar.gz

RUN ln -s ncbi-blast-${BLAST_VERSION}+ blast

ENV PATH "$PATH:/opt/blast/bin:/opt/blast"


CMD ["/bin/bash"]