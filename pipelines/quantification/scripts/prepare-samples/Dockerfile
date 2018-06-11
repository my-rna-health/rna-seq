FROM quay.io/comp-bio-aging/pipeline-script:latest

LABEL maintainer="antonkulaga@gmail.com"

ENV JAVA_OPTS="-Xmx3072m"

ADD scripts /scripts/
RUN chmod +x /scripts/*
WORKDIR /scripts

RUN /scripts/tsv.sc info
RUN /scripts/ml.sc info
RUN /scripts/process.sc info

CMD /scripts/process.sc