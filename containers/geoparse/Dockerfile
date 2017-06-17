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

WORKDIR /data

ENTRYPOINT ["/opt/geoparse/run.py"]

CMD ["--help"]
