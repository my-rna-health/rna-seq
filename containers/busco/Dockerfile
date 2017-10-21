FROM quay.io/comp-bio-aging/hmmer:latest

MAINTAINER Anton Kulaga

USER root

RUN apt-get update && apt-get upgrade -y -q
RUN apt-get install -y libsqlite3-dev libboost-graph-dev liblpsolve55-dev curl

# install ncbi blast
RUN cd /root && wget -O- https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/ncbi-blast-2.6.0+-x64-linux.tar.gz | tar zx  && \
 cp ncbi-blast*/bin/* /usr/bin/

RUN rm /usr/bin/python
RUN ln -s /usr/bin/python3 /usr/bin/python

WORKDIR /opt

## -- install busco -- ##
RUN git clone --recursive https://gitlab.com/ezlab/busco.git
WORKDIR /opt/busco
RUN curl -SL http://busco.ezlab.org/datasets/vertebrata_odb9.tar.gz | tar -xz \
	&& curl -SL http://busco.ezlab.org/datasets/metazoa_odb9.tar.gz | tar -xz \
	&& curl -SL http://busco.ezlab.org/datasets/eukaryota_odb9.tar.gz | tar -xz \
	&& curl -SL http://busco.ezlab.org/datasets/arthropoda_odb9.tar.gz | tar -xz \
	&& curl -SL http://busco.ezlab.org/datasets/fungi_odb9.tar.gz | tar -xz \
	&& curl -SL http://busco.ezlab.org/datasets/bacteria_odb9.tar.gz | tar -xz \
	&& python3 setup.py install

# Use the default config file
RUN mv /opt/busco/config/config.ini.default /opt/busco/config/config.ini

ENV PATH "$PATH:/opt/busco:/opt/busco/scripts"
CMD ["run_BUSCO.py"]
