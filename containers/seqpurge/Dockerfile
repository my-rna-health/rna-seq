FROM ubuntu:latest
MAINTAINER Anton Kulaga

RUN apt-get update && \
    apt-get install -y git build-essential g++ qt5-default libqt5xmlpatterns5-dev libqt5sql5-mysql git cmake python python-matplotlib

RUN cd /root/ && git clone --recursive https://github.com/imgag/ngs-bits && cd ngs-bits && \
    sed -i "s/[^qc]make;/ make -j$(nproc);/g" Makefile && sed -i "s/[^qc]make$/ make -j$(nproc)/g" Makefile && \
    make build_3rdparty && make build_tools_release

WORKDIR /root/ngs-bits

RUN cp bin/SeqPurge /usr/bin/ && \
    cp bin/libcppCORE.so.1 bin/libcppNGS.so.1 bin/libcppXML.so.1 bin/libbamtools.so.2.4.1 /usr/lib/ && \
    rm -rf /root/ngs-bits

RUN mkdir /data
WORKDIR /data
CMD ["SeqPurge"]