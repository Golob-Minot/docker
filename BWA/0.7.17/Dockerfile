
# # bwa
#
# VERSION               golob/bwa:0.7.17__bcw.0.3.0I

FROM      ubuntu:16.04
RUN export DEBIAN_FRONTEND=noninteractive

RUN apt-get update && apt-get install -y \
python3-dev \
python3-pip \
wget \
perl \
zlib1g-dev \
libbz2-dev \
liblzma-dev \
libcurl3-dev \
libncursesw5-dev \
pigz \
awscli \
&& apt-get clean \
&& apt-get purge \
&& rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/*

RUN ln -s /usr/bin/python3 /usr/local/bin/python
RUN pip3 install \
bucket_command_wrapper==0.3.0

RUN mkdir /src
WORKDIR /src
RUN wget https://sourceforge.net/projects/bio-bwa/files/bwa-0.7.17.tar.bz2/download && \
tar xjvf download && rm download
WORKDIR /src/bwa-0.7.17
RUN make
RUN cp bwa /usr/local/bin/ && \
cp qualfa2fq.pl /usr/local/bin && \
cp xa2multi.pl /usr/local/bin

WORKDIR /src
RUN wget https://sourceforge.net/projects/samtools/files/samtools/1.9/samtools-1.9.tar.bz2/download
RUN tar xjvf download
WORKDIR /src/samtools-1.9
RUN ./configure
RUN make && make install

WORKDIR /
RUN rm -r /src/

RUN mkdir -p ~/.aws
RUN mkdir -p /mnt/inputs/s3 && mkdir -p /mnt/inputs/file 
RUN mkdir -p /mnt/outputs/s3 && mkdir -p /mnt/outputs/file 
RUN mkdir /scratch && mkdir /working
RUN mkdir /fh/ && mkdir /app