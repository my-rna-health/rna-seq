FROM ubuntu:latest

MAINTAINER Anton Kulaga <antonkulaga@gmail.com>

RUN apt-get update && \
 apt-get upgrade -y && \
 apt-get install -y cmake zip wget gcc git build-essential

RUN apt-get install -y \
    libboost-iostreams-dev zlib1g-dev libsqlite3-dev \
    libgsl-dev libsuitesparse-dev liblpsolve55-dev \
    bamtools libbamtools-dev

WORKDIR /opt

RUN wget http://bioinf.uni-greifswald.de/augustus/binaries/augustus.current.tar.gz && \
    tar -xzvf augustus* && \
    rm -rf *.tar.gz

WORKDIR /opt/augustus
RUN make

ENV AUGUSTUS_CONFIG_PATH /opt/augustus/config
ENV PATH /usr/local/sbin:/usr/local/bin:/usr/sbin:/usr/bin:/sbin:/bin:/opt/augustus/bin:/opt/augustus/scripts
CMD ["/opt/augustus/bin/augustus"]