FROM ubuntu:latest

RUN sudo apt-get install zlib1g-dev libbz2-dev
RUN sudo apt-get install build-essential

RUN wget https://github.com/MikkelSchubert/adapterremoval/archive/v2.2.2.tar.gz -O adapterremoval-2.2.2.tar.gz \
    && tar zxf adapterremoval-2.2.2.tar.gz  \
    && mkdir /opt/adapterremoval \
    && ln -s /adapterremoval-2.2.2 /opt/adapterremoval
WORKDIR /opt/adapterremoval
RUN make

RUN mkdir /data
WORKDIR /opt

CMD ["--help"]
