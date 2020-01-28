# barcodecop
#
# VERSION               golob/barcodecop:0.5__bc

FROM      biocontainers/biocontainers:vdebian-buster-backports_cv1
ENV DEBIAN_FRONTEND noninteractive
USER root
RUN apt-get update && apt-get install -y wget python3-pip \ 
&& pip3 install barcodecop \
&& apt-get clean && apt-get purge && rm -rf /var/lib/apt/lists/* /tmp/*
