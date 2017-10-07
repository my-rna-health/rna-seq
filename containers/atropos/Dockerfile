FROM snakeego/cython

MAINTAINER Anton Kulaga <antonkulaga@gmail.com>

# bash support to enable container usage from cromwell and other workflow engines
RUN apk add --update bash && rm -rf /var/cache/apk/*

# GCC
RUN export NPROC=$(grep -c ^processor /proc/cpuinfo 2>/dev/null || 1) \
    && apk --no-cache add --virtual build-deps \
        musl-dev \
        linux-headers \
        g++ \
        zlib-dev bzip2-dev xz-dev \
        git

# Additional Atropos dependencies
RUN pip install tqdm pysam jinja2

# git clone atropos and install
WORKDIR /opt
RUN git clone https://github.com/jdidion/atropos.git
RUN cd /opt/atropos/ && make install

# Cleanup
RUN rm -rf /var/cache/apk/* /atropos

ENV PATH "$PATH:/opt/atropos"

RUN mkdir /data
WORKDIR /data

CMD ["atropos"]