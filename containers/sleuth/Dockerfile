FROM bioconductor/release_base2

MAINTAINER maintainer@bioconductor.org

ADD install.R /tmp/

# invalidates cache every 24 hours
ADD http://master.bioconductor.org/todays-date /tmp/

RUN R -f /tmp/install.R