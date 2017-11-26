FROM quay.io/comp-bio-aging/blast:latest

MAINTAINER Anton Kulaga

WORKDIR /opt

## -- install busco -- ##
RUN git clone --recursive https://gitlab.com/ezlab/busco.git
WORKDIR /opt/busco
RUN curl -SL http://busco.ezlab.org/datasets/vertebrata_odb9.tar.gz | tar -xz \
	&& curl -SL http://busco.ezlab.org/datasets/metazoa_odb9.tar.gz | tar -xz \
	&& curl -SL http://busco.ezlab.org/datasets/eukaryota_odb9.tar.gz | tar -xz \
	&& curl -SL http://busco.ezlab.org/datasets/arthropoda_odb9.tar.gz | tar -xz \
	&& curl -SL http://busco.ezlab.org/datasets/fungi_odb9.tar.gz | tar -xz \
	&& curl -SL http://busco.ezlab.org/datasets/bacteria_odb9.tar.gz | tar -xz

RUN /usr/bin/python3 setup.py install

# Use the default config file
ADD config.ini /opt/busco/config/

ENV PATH "$PATH:/opt/busco:/opt/busco/scripts"

#ADDING R FOR SCRIPTING
RUN apt-get update && apt-get -y install r-base

#setup R configs
RUN echo "r <- getOption('repos'); r['CRAN'] <- 'http://cran.us.r-project.org'; options(repos = r);" > ~/.Rprofile
RUN Rscript -e "install.packages('ggplot2')"

CMD ["run_BUSCO.py"]