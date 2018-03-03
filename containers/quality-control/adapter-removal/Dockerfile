FROM phusion/baseimage:latest
WORKDIR /tmp
ENV VERSION v2.2.2
ENV BUILD_PACKAGES \
    build-essential \
    g++ \
    libbz2-dev \
    zlib1g-dev \
    git \
    perl
ENV RUNTIME_PACKAGES \
    time
RUN apt-get update \
    && DEBIAN_FRONTEND=noninteractive apt-get -y upgrade \
    && DEBIAN_FRONTEND=noninteractive apt-get install -y -qq \
        $BUILD_PACKAGES \
        $RUNTIME_PACKAGES \
    && git clone --recursive --branch $VERSION https://github.com/MikkelSchubert/adapterremoval \
    && cd adapterremoval \
    && make \
    && mv /tmp/adapterremoval/build/AdapterRemoval /usr/local/bin \
    && cd .. \
    && apt-get remove --purge -y $BUILD_PACKAGES \
    && apt-get autoremove -y \
    && rm -rf /var/lib/apt/lists/* /tmp/* \
    && locale-gen en_US.UTF-8 \
    && update-locale LANG=en_US.UTF-8


CMD ["AdapterRemoval"]
