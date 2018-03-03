FROM quay.io/comp-bio-aging/base:latest

MAINTAINER Anton Kulaga

ENV DENOTATE="1.11"

WORKDIR /opt/

RUN wget http://deweylab.biostat.wisc.edu/detonate/detonate-${DENOTATE}-precompiled.tar.gz && \
    tar -zxvf detonate-${DENOTATE}-precompiled.tar.gz && \
    rm -rf detonate-${DENOTATE}-precompiled.tar.gz

RUN ln -s detonate-${DENOTATE}-precompiled denotate

WORKDIR /opt/denotate

ENV PATH "$PATH:/opt/denotate:/opt/denotate"

RUN mkdir /data
WORKDIR /data

CMD ["/bin/bash"]