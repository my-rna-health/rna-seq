FROM quay.io/comp-bio-aging/pipeline-script:latest

LABEL maintainer="antonkulaga@gmail.com"

ADD scripts /scripts/
RUN chmod +x /scripts/*
RUN /usr/local/bin/amm /scripts/dependencies.sc

CMD /scripts/run.sc