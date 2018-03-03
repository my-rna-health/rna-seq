FROM ruby:2.2.8

MAINTAINER Anton Kulaga <antonkulaga@gmail.com>

RUN apt-get update && \
 apt-get upgrade -y && \
 apt-get install -y cmake zip wget gcc git build-essential

ENV TRANSRATE_VERSION="1.0.3"
ENV TRANSRATE_PLATFORM="linux-x86_64"

WORKDIR /opt/

#RUN wget https://bintray.com/blahah/generic/download_file?file_path=transrate-${TRANSRATE_VERSION}-${TRANSRATE_PLATFORM}.tar.gz && \
#    tar -zxvf transrate-${TRANSRATE_VERSION}-${TRANSRATE_PLATFORM}.tar.gz && \
#    rm -rf transrate-${TRANSRATE_VERSION}-${TRANSRATE_PLATFORM}.tar.gz && \
#    ln -s transrate-${TRANSRATE_VERSION}-${TRANSRATE_PLATFORM} transrate

RUN gem install transrate
RUN transrate --install-deps all

WORKDIR /opt/transrate

RUN mkdir /data
WORKDIR /data

CMD ["transrate"]

