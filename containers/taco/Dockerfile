FROM ubuntu:latest
MAINTAINER Anton Kulaga

ENV VERSION="0.7.3"
ENV PLATFORM="Linux_x86_64"

RUN apt-get update && \
 apt-get upgrade -y && \
 apt-get install -y cmake zip wget gcc git build-essential

WORKDIR /opt/

RUN wget https://github.com/tacorna/taco/releases/download/v$VERSION/taco-v$VERSION.$PLATFORM.tar.gz && \
    tar -zxvf taco-v$VERSION.$PLATFORM.tar.gz && \
    rm -rf taco-v$VERSION.$PLATFORM.tar.gz && \
    ln -s taco-v$VERSION.$PLATFORM taco

WORKDIR /opt/taco

ENV PATH "$PATH:/opt/taco"

RUN mkdir /data
WORKDIR /data

CMD ["/bin/bash"]