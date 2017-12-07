FROM quay.io/comp-bio-aging/base

MAINTAINER Anton Kulaga <antonkulaga@gmail.com>

ENV DIAMOND_VERSION="v0.9.13"
ENV DIAMOND_PLATFORM="linux64"

ENV DIAMOND_PATH /opt/diamond-$DIAMOND_PLATFORM
RUN mkdir $DIAMOND_PATH

WORKDIR $DIAMOND_PATH

RUN wget http://github.com/bbuchfink/diamond/releases/download/$DIAMOND_VERSION/diamond-$DIAMOND_PLATFORM.tar.gz && \
    tar -xzvf diamond-$DIAMOND_PLATFORM.tar.gz && \
    rm -rf *.tar.gz

ENV DIAMOND_PATH /opt/diamond-$DIAMOND_PLATFORM

ENV PATH $PATH:$DIAMOND_PATH
CMD ["diamond"]