FROM ubuntu:17.10

MAINTAINER Anton Kulaga <antonkulaga@gmail.com>

RUN apt-get update && \
    apt-get install -y \
      libncurses5-dev \
      libncursesw5-dev \
      build-essential \
      zlib1g-dev \
      libbz2-dev \
      liblzma-dev \
      wget && \
    apt-get clean && \
    apt-get purge && \
    rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/*
    
    
    
ENV VERSION="1.6"
ENV PLATFORM="linux-intel-x86_64"

WORKDIR /opt/

RUN wget https://github.com/samtools/samtools/releases/download/$VERSION/samtools-$VERSION.tar.bz2 && \
    tar -xjf samtools-$VERSION.tar.bz2 && \
    rm -rf samtools-$VERSION.tar.bz2 && \
    ln -s samtools-$VERSION samtools

WORKDIR /opt/samtools

RUN ./configure
RUN make


ENV PATH "$PATH:/opt/samtools"

CMD ["/bin/bash"]