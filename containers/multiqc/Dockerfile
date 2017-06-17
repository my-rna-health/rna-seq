FROM snakeego/cython

MAINTAINER Anton Kulaga antonkulaga@gmail.com

# bash support to enable container usage from cromwell and other workflow engines
RUN apk add --update bash && rm -rf /var/cache/apk/*

# GCC
RUN export NPROC=$(grep -c ^processor /proc/cpuinfo 2>/dev/null || 1) \
    && apk --no-cache add --virtual build-deps \
        musl-dev \
        linux-headers \
        g++ \
        zlib-dev bzip2-dev xz-dev

RUN python -m pip install --upgrade pip

RUN pip install multiqc click

RUN mkdir /opt
RUN mkdir /opt/multiqc
RUN mkdir /data

COPY src/ /opt/multiqc

WORKDIR /opt/geoparse
RUN python setup.py build_ext --inplace
RUN chmod +x /opt/geoparse/run.py

ENV PATH "$PATH:/opt/multiqc"

WORKDIR /data

ENTRYPOINT ["/opt/multiqc/run.py"]

CMD ["--help"]