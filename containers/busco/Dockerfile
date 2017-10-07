FROM blaxterlab/braker

MAINTAINER Anton Kulaga

ENV AUGUSTUS_CONFIG_PATH /opt/augustus-3.2.2/config

RUN apt-get update && \
 apt-get upgrade -y && \
 apt-get install -y wget zlib1g-dev libboost-iostreams-dev libsqlite3-dev libboost-graph-dev liblpsolve55-dev cmake gcc git build-essential

# install hmmer
RUN cd  && wget -O- http://eddylab.org/software/hmmer3/3.1b2/hmmer-3.1b2-linux-intel-x86_64.tar.gz | tar zx && \
 cd hmmer-3.1b2-linux-intel-x86_64/ && ./configure && make && make install

# install ncbi blast
RUN cd /root && wget -O- https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/ncbi-blast-2.6.0+-x64-linux.tar.gz | tar zx  && \
 cp ncbi-blast*/bin/* /usr/bin/

# install busco
RUN cd /root && git clone http://gitlab.com/ezlab/busco && cp busco/*.py /usr/bin/

# make python3 the default python
RUN ln -s /usr/bin/python3 /usr/bin/python

WORKDIR /data

ENTRYPOINT ["BUSCO.py"]

CMD ["--help"]