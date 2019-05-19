FROM quay.io/comp-bio-aging/ammonite:latest

LABEL maintainer="antonkulaga@gmail.com"

ENV JAVA_OPTS="-Xmx3072m"

ADD scripts /scripts/
RUN chmod +x /scripts/*
WORKDIR /scripts
RUN amm dependencies.sc

CMD /scripts/sample_table.sc