## golob/ngmerge:0.3__bcw.0.3.1

FROM      ubuntu:18.04
# For singularity on the hutch cluster
ENV DEBIAN_FRONTEND=noninteractive

RUN apt-get update && apt-get install -y \
python-dev \
python-pip \
wget \
pigz \
awscli \
zlib1g-dev \
&& apt-get clean \
&& apt-get purge \
&& rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/*

RUN mkdir -p /src/ && \
cd /src/ && \
wget https://github.com/harvardinformatics/NGmerge/archive/v0.3.tar.gz && \
tar xzvf v0.3.tar.gz && \
cd /src/NGmerge-0.3/ && \
make && \
cp NGmerge /usr/local/bin/ && \
cd / && \
rm -rf /src/

RUN pip install bucket_command_wrapper==0.3.1

RUN mkdir -p /mnt/inputs/file && mkdir -p /mnt/outputs/file
# Various directories
RUN mkdir -p ~/.aws
RUN mkdir -p /mnt/inputs/s3 && mkdir -p /mnt/inputs/file 
RUN mkdir -p /mnt/outputs/s3 && mkdir -p /mnt/outputs/file 
RUN mkdir /scratch && mkdir /working
RUN mkdir /fh/ && mkdir /app