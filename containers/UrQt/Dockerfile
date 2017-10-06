FROM ubuntu:latest

MAINTAINER Anton Kulaga <antonkulaga@gmail.com>

ENV VERSION 1.0.17

RUN apt-get update && apt-get install -qqy build-essential zlib1g-dev wget

WORKDIR /opt
RUN wget ftp://pbil.univ-lyon1.fr/pub/logiciel/UrQt/UrQt.$VERSION.tar.gz
RUN tar xvzf UrQt.$VERSION.tar.gz
RUN ln -s /opt/UrQt.$VERSION /opt/UrQt

WORKDIR /opt/UrQt
RUN make


ENV PATH "$PATH:/opt/UrQt"

RUN mkdir /data
WORKDIR /data

CMD ["/opt/UrQt/UrQt"]
