FROM ubuntu:17.10

MAINTAINER Anton Kulaga <antonkulaga@gmail.com>

RUN apt-get update && \
    apt-get upgrade -y && \
    apt-get install -y python3-software-properties cmake zip wget gcc git build-essential

RUN apt-get install -y \
    libboost-iostreams-dev zlib1g-dev libsqlite3-dev \
    bamtools libbamtools-dev python3



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