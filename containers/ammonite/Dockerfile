FROM quay.io/comp-bio-aging/java:latest

LABEL maintainer "antonkulaga@gmail.com"

RUN curl -L -o /usr/local/bin/amm https://git.io/vASZm && chmod +x /usr/local/bin/amm && amm

RUN mkdir /data

CMD ["/usr/local/bin/amm"]