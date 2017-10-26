FROM ubuntu:17.10

MAINTAINER Anton Kulaga <antonkulaga@gmail.com>

RUN apt-get update && \
    apt-get upgrade -y && \
    apt-get install -y python3-software-properties cmake zip wget gcc git build-essential curl

RUN ln -s /usr/bin/python3 /usr/bin/python

# Install Java.
RUN \
  echo oracle-java8-installer shared/accepted-oracle-license-v1-1 select true | debconf-set-selections && \
  add-apt-repository -y ppa:webupd8team/java && \
  apt-get update && \
  apt-get install -y oracle-java8-installer && \
  rm -rf /var/lib/apt/lists/* && \
  rm -rf /var/cache/oracle-jdk8-installer


# Define working directory.
WORKDIR /data

# Define commonly used JAVA_HOME variable
ENV JAVA_HOME /usr/lib/jvm/java-8-oracle

# Define default command.
CMD ["bash"]