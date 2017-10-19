FROM python:3.6.3

LABEL \
  author="Anton Kulaga" \
  description="MultiQC" \
  maintainer="antonkulaga@gmail.com"

# Install libraries
RUN \
  apt-get update && apt-get install -y --no-install-recommends \
  g++ \
  git \
  wget \
  && wget --quiet -O /opt/get-pip.py https://bootstrap.pypa.io/get-pip.py \
  && python /opt/get-pip.py \
  && rm -rf /var/lib/apt/lists/* /opt/get-pip.py

# Install MultiQC
RUN pip install git+git://github.com/ewels/MultiQC.git
CMD ["multiqc"]