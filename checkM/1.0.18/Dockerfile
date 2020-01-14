# checkM
#
# VERSION               golob/checkm:1.0.18

FROM      ubuntu:18.04

ENV DEBIAN_FRONTEND=noninteractive
RUN apt-get update && apt-get install -y \
python-dev \
python-pip \
zlib1g-dev \
libbz2-dev \
liblzma-dev \
build-essential \
wget \
unzip \
&& apt-get clean \
&& apt-get purge \
&& rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/*

# Install pplacer
RUN mkdir /src/
COPY pplacer-linux-v1.1.alpha19-4-g1189285.zip /src/
RUN cd /src/ && \
unzip pplacer-linux-v1.1.alpha19-4-g1189285.zip && \
cp guppy /usr/local/bin && \
cp pplacer /usr/local/bin && \
cp rppr /usr/local/bin && \
cd /src/scripts && python setup.py install && \
rm -r /src/*

# Install prodigal 2.6.3
RUN cd /usr/local/bin/ && \
wget https://github.com/hyattpd/Prodigal/releases/download/v2.6.3/prodigal.linux && \
mv prodigal.linux prodigal && \
chmod +x prodigal

# Install HMMER 3.2.1
RUN cd /src && \
wget http://eddylab.org/software/hmmer/hmmer-3.2.1.tar.gz && \
tar xzvf hmmer-3.2.1.tar.gz && \
cd /src/hmmer-3.2.1 && \
./configure && \
make && make check && make install && \
rm -r /src/*

RUN pip install \
awscli \
numpy \
pysam \
biopython \
checkm-genome

# Get the checkM database
RUN mkdir -p /db/checkm && \
cd /db/checkm && \
wget https://data.ace.uq.edu.au/public/CheckM_databases/checkm_data_2015_01_16.tar.gz && \
tar xzvf checkm_data_2015_01_16.tar.gz && \
rm checkm_data_2015_01_16.tar.gz && \
checkm data setRoot /db/checkm

RUN mkdir /scratch
WORKDIR /scratch

