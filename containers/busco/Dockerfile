FROM blaxterlab/braker

MAINTAINER Anton Kulaga

USER root

RUN apt-get update && apt-get upgrade -y -q
RUN apt-get install -y libsqlite3-dev libboost-graph-dev liblpsolve55-dev curl

# install hmmer
RUN cd  && wget -O- http://eddylab.org/software/hmmer3/3.1b2/hmmer-3.1b2-linux-intel-x86_64.tar.gz | tar zx && \
 cd hmmer-3.1b2-linux-intel-x86_64/ && ./configure && make && make install

# install ncbi blast
RUN cd /root && wget -O- https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/ncbi-blast-2.6.0+-x64-linux.tar.gz | tar zx  && \
 cp ncbi-blast*/bin/* /usr/bin/

WORKDIR /opt

## -- install busco -- ##
RUN git clone --recursive https://gitlab.com/ezlab/busco.git
WORKDIR /opt/busco
RUN curl -SL http://busco.ezlab.org/v2/datasets/vertebrata_odb9.tar.gz | tar -xz \
	&& curl -SL http://busco.ezlab.org/v2/datasets/metazoa_odb9.tar.gz | tar -xz \
	&& curl -SL http://busco.ezlab.org/v2/datasets/eukaryota_odb9.tar.gz | tar -xz \
	&& curl -SL http://busco.ezlab.org/v2/datasets/arthropoda_odb9.tar.gz | tar -xz \
	&& curl -SL http://busco.ezlab.org/v2/datasets/fungi_odb9.tar.gz | tar -xz \
	&& curl -SL http://busco.ezlab.org/v2/datasets/bacteria_odb9.tar.gz | tar -xz \
	&& python setup.py install

RUN rm /usr/bin/python
RUN ln -s /usr/bin/python3 /usr/bin/python

ENV PATH "$PATH:/opt/busco"
CMD ["BUSCO.py"]
