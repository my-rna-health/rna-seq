FROM quay.io/comp-bio-aging/scala:latest

LABEL maintainer "antonkulaga@gmail.com"

RUN curl -L -o /usr/local/bin/amm https://git.io/vASZm && chmod +x /usr/local/bin/amm

RUN echo 'println("Hello ammonite!")' > hw.scala /usr/local/bin/amm hw.scala 2>/dev/null a && rm -rf hw.scala

CMD ["amm"]