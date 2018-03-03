FROM quay.io/comp-bio-aging/java:latest

LABEL maintainer "antonkulaga@gmail.com"

# Env variables
ENV SBT_VERSION 1.1.1

# Install sbt
RUN \
  curl -L -o sbt-$SBT_VERSION.deb https://dl.bintray.com/sbt/debian/sbt-$SBT_VERSION.deb && \
  dpkg -i sbt-$SBT_VERSION.deb && \
  rm sbt-$SBT_VERSION.deb && \
  apt-get update && \
  apt-get install sbt && \
  sbt sbtVersion

CMD [sbt]