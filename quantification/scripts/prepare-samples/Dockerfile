FROM quay.io/comp-bio-aging/pipeline-script:latest

LABEL maintainer="antonkulaga@gmail.com"

ADD scripts /scripts/
RUN chmod +x /scripts/*
RUN /scripts/run.sc info

CMD /scripts/run.sc