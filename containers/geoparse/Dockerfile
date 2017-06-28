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

ENV PATH "$PATH:/opt/geoparse"

# sratoolkit installation
WORKDIR /opt
RUN \
	wget -q https://ftp-trace.ncbi.nlm.nih.gov/sra/sdk/2.8.2/sratoolkit.2.8.2-ubuntu64.tar.gz -O sratoolkit.2.8.2.tar.gz \
    && tar zxf sratoolkit.2.8.2.tar.gz  \
    && mkdir /opt/sratoolkit \
    && ln -s /opt/sratoolkit.2.8.2-ubuntu64/bin/* /opt/sratoolkit/

# set environment
ENV PATH="$PATH:/opt/sratoolkit/"

WORKDIR /data

ENTRYPOINT ["/opt/geoparse/run.py"]

CMD ["--help"]
