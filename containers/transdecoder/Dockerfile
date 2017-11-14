FROM bioperl/bioperl

MAINTAINER Anton Kulaga antonkulaga@gmail.com

ENV VERSION 5.0.2

RUN apt-get update && \
 apt-get upgrade -y && \
 apt-get install -y cmake zip wget gcc git build-essential

WORKDIR /opt/

RUN wget https://github.com/TransDecoder/TransDecoder/archive/v$VERSION.zip && \
    unzip v$VERSION.zip && \
    ln -s TransDecoder-$VERSION TransDecoder

WORKDIR /opt/TransDecoder

ENV PATH "$PATH:/opt/TransDecoder"

RUN mkdir /data
WORKDIR /data

CMD ["/bin/bash"]
