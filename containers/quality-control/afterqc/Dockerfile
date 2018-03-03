FROM ubuntu:latest

MAINTAINER Anton Kulaga <antonkulaga@gmail.com>

RUN apt-get update && apt-get install -y ca-certificates software-properties-common

# install python
RUN apt-get install -y python-pip python python-dev wget

# install pypy
RUN add-apt-repository ppa:pypy/ppa
RUN apt-get update
RUN apt-get install -y pypy pypy-dev

ENV VERSION 0.9.6

WORKDIR /opt
RUN wget https://github.com/OpenGene/AfterQC/archive/v$VERSION.tar.gz
RUN tar xvzf v$VERSION.tar.gz
RUN ln -s /opt/AfterQC-$VERSION /opt/AfterQC

WORKDIR /opt/AfterQC
RUN make

ENV PATH "$PATH:/opt/AfterQC"

CMD ["/opt/AfterQC/after.py"]
