FROM python:3

MAINTAINER Anton Kulaga antonkulaga@gmail.com

#installing mypy
RUN python3 -m pip install -U mypy

COPY requirements.txt ./
RUN pip install --no-cache-dir -r requirements.txt

RUN mkdir /opt/geoparse
RUN mkdir /data

COPY src/ /opt/geoparse

WORKDIR /opt/geoparse
RUN chmod +x /opt/geoparse/run.py

ENV PLATFORM="ubuntu64"
ENV SRA_VERSION="2.8.2-1"
#ENV SRA_VERSION="2.9.0"

ENV PATH "$PATH:/opt/geoparse"

# sratoolkit installation
WORKDIR /opt
RUN \
	wget -q https://ftp-trace.ncbi.nlm.nih.gov/sra/sdk/${SRA_VERSION}/sratoolkit.${SRA_VERSION}-${PLATFORM}.tar.gz -O sratoolkit.${SRA_VERSION}.tar.gz \
    && tar zxf sratoolkit.${SRA_VERSION}.tar.gz  \
    && mkdir /opt/sratoolkit \
    && ln -s /opt/sratoolkit.${SRA_VERSION}-${PLATFORM}/bin/* /opt/sratoolkit/

# set environment
ENV PATH="$PATH:/opt/sratoolkit/"

WORKDIR /data

CMD ["/opt/geoparse/run.py"]