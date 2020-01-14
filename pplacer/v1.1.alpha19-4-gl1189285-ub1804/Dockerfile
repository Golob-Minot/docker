# pplacer
#
# VERSION               golob/pplacer:1.1alpha19rc_BCW_0.3.1A

FROM      ubuntu:18.04

ENV DEBIAN_FRONTEND=noninteractive
RUN apt-get update && apt-get install -y \
python-dev \
python3-dev \
python3-pip \
postgresql \
libpq-dev \
wget \
unzip \
&& apt-get clean \
&& apt-get purge \
&& rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/*

RUN pip3 install \
numpy \
biopython \
psycopg2 \
awscli \
boto3 \
bucket_command_wrapper \
taxtastic \
csvkit

RUN mkdir /fh && mkdir /app && mkdir /src && mkdir /scratch
COPY pplacer-linux-v1.1.alpha19-4-g1189285.zip /src/
RUN cd /src/ && \
unzip pplacer-linux-v1.1.alpha19-4-g1189285.zip && \
cp guppy /usr/local/bin && \
cp pplacer /usr/local/bin && \
cp rppr /usr/local/bin && \
cd /src/scripts && python setup.py install && \
rm -r /src/

RUN mkdir -p /pplacer/refpkg && chmod a+w /pplacer/refpkg
WORKDIR /scratch

