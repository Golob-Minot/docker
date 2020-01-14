# spades
#
# VERSION               golob/spades:3.12.0__bcw.0.3.0B

FROM      ubuntu:16.04
# For singularity on the hutch cluster
RUN export DEBIAN_FRONTEND=noninteractive
RUN mkdir /fh && mkdir /app && mkdir /src
RUN mkdir -p /mnt/inputs/file && mkdir -p /mnt/outputs/file && mkdir /scratch && mkdir /working
RUN apt-get update && apt-get install -y \
python3-dev \
python3-pip \
wget \
cmake \
zlib1g-dev \
libbz2-dev

RUN ln -s /usr/bin/python3 /usr/local/bin/python
RUN pip3 install \
awscli>=1.15.14 \
boto3>=1.7.14 \
bucket_command_wrapper==0.3.0 

# Get and extract spades
WORKDIR /src/
#RUN wget http://cab.spbu.ru/files/release3.12.0/SPAdes-3.12.0.tar.gz
#RUN tar xzvf SPAdes-3.12.0.tar.gz
#WORKDIR /src/SPAdes-3.12.0
#RUN PREFIX=/usr/local ./spades_compile.sh
#WORKDIR /scratch
#RUN spades.py --test
#RUN rm -r /src/
##
RUN wget http://cab.spbu.ru/files/release3.12.0/SPAdes-3.12.0-Linux.tar.gz && \
tar xzvf SPAdes-3.12.0-Linux.tar.gz && \
cp /src/SPAdes-3.12.0-Linux/bin/* /usr/local/bin/ && \
cp -a /src/SPAdes-3.12.0-Linux/share/spades /usr/local/share/

WORKDIR /


