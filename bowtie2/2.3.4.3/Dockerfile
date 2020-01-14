# # bowtie2
#
# VERSION               golob/bowtie2:2.3.4.3__bcw.0.3.0

FROM      ubuntu:18.04
RUN export DEBIAN_FRONTEND=noninteractive
RUN mkdir /fh && mkdir /app && mkdir /src
RUN mkdir -p /mnt/inputs/file && mkdir -p /mnt/outputs/file && mkdir /scratch && mkdir /working
RUN apt-get update && apt-get install -y \
python3-dev \
python3-pip \
wget \
zlib1g-dev \
libbz2-dev
RUN ln -s /usr/bin/python3 /usr/local/bin/python

WORKDIR /src/
RUN wget https://sourceforge.net/projects/bowtie-bio/files/bowtie2/2.3.4.3/bowtie2-2.3.4.3-linux-x86_64.zip/download
RUN unzip download
RUN cp /src/bowtie2-2.3.4.3-linux-x86_64/bowtie* /usr/local/bin/
WORKDIR /