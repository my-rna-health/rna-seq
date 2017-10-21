FROM ruby:2.2.8

MAINTAINER Anton Kulaga <antonkulaga@gmail.com>

RUN apt-get update && \
 apt-get upgrade -y && \
 apt-get install -y cmake zip wget gcc git build-essential

RUN gem install transrate
RUN transrate --install-deps all

RUN mkdir /data
WORKDIR /data

CMD ["transrate"]

