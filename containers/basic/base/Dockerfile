FROM ubuntu:17.10

MAINTAINER Anton Kulaga <antonkulaga@gmail.com>

RUN apt-get update && \
    apt-get upgrade -y && \
    apt-get install -y python3-software-properties cmake zip wget gcc git build-essential curl

RUN ln -s /usr/bin/python3 /usr/bin/python

# Define default command.
CMD ["/bin/bash"]