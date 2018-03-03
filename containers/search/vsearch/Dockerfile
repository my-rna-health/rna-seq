FROM quay.io/comp-bio-aging/base:latest

MAINTAINER Anton Kulaga

ENV PLATFORM=linux-x86_64

ENV VSEARCH_VERSION=2.6.2

WORKDIR /opt

RUN wget https://github.com/torognes/vsearch/releases/download/v${VSEARCH_VERSION}/vsearch-${VSEARCH_VERSION}-${PLATFORM}.tar.gz && \
     tar -zxvf vsearch-${VSEARCH_VERSION}-${PLATFORM}.tar.gz && \
     rm -rf vsearch-${VSEARCH_VERSION}-${PLATFORM}.tar.gz

RUN ln -s vsearch-${VSEARCH_VERSION}-${PLATFORM} vsearch

ENV PATH "$PATH:/opt/vsearch/bin"

CMD ["/opt/vsearch/bin/vsearch", "--help"]