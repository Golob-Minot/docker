# spades
#
# VERSION               golob/spades:3.13.1__bcw.0.3.1

FROM      ubuntu:16.04
# For singularity on the hutch cluster
ENV DEBIAN_FRONTEND=noninteractive
RUN apt-get update && apt-get install -y \
python3-dev \
python3-pip \
wget \
cmake \
zlib1g-dev \
libbz2-dev \
&& apt-get clean \
&& apt-get purge \
&& rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/*

RUN ln -s /usr/bin/python3 /usr/local/bin/python
RUN pip3 install \
awscli>=1.15.14 \
boto3>=1.7.14 \
bucket_command_wrapper==0.3.1

# Get and extract spades

RUN mkdir /src && cd /src/ && \
wget https://github.com/ablab/spades/releases/download/v3.13.1/SPAdes-3.13.1-Linux.tar.gz && \
tar xzvf SPAdes-3.13.1-Linux.tar.gz && \
cp /src/SPAdes-3.13.1-Linux/bin/* /usr/local/bin/ && \
cp -a /src/SPAdes-3.13.1-Linux/share/spades /usr/local/share/ && \
cd / && rm -r /src/

RUN mkdir /fh && mkdir /app
RUN mkdir -p /mnt/inputs/file && mkdir -p /mnt/outputs/file && mkdir /scratch && mkdir /working

WORKDIR /


