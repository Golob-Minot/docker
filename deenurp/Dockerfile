# golob/deenurp:0.2.6
FROM      ubuntu:18.04

env DEBIAN_FRONTEND=noninteractive
RUN apt-get update && apt-get install -y \
python-dev \
python-pip \
wget \
&& apt-get clean \
&& apt-get purge \
&& rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/*

RUN pip install \
setuptools==38.4.0 \
numpy==1.14.0 \
biopython==1.66

RUN pip install \
MarkupSafe==1.0 \
six==1.11.0 \ 
DendroPy==4.3.0 \
Jinja2==2.10 \
PyYAML==3.12 \
SQLAlchemy==1.2.0 \
decorator==4.1.2 \
fastalite==0.3 \
psycopg2==2.7.3.2 \
python-dateutil==2.6.1 \
pytz==2017.3 \ 
scikit-learn==0.19.1 \
Cython==0.27.3 \
futures==3.2.0 \
hdbscan==0.8.11 \
pandas==0.22.0 \
peasel==0.3.1 \
scipy==1.0.0 \ 
seqmagick==0.6.2 \
taxtastic==0.8.5

RUN pip install awscli

RUN mkdir -p /src/ && cd /src && \
wget https://github.com/fhcrc/deenurp/archive/v0.2.6.tar.gz && \
tar xzvf v0.2.6.tar.gz && \
cd /src/deenurp-0.2.6/ && python setup.py install && \
rm -r /src/*

RUN mkdir -p /src/ && cd /src && \
wget http://eddylab.org/infernal/infernal-1.1.2-linux-intel-gcc.tar.gz && \
tar xzvf infernal-1.1.2-linux-intel-gcc.tar.gz && \
cp /src/infernal-1.1.2-linux-intel-gcc/binaries/* /usr/local/bin/ && \
rm -r /src/*

RUN mkdir -p /src/ && cd /usr/local/bin && \
wget http://www.microbesonline.org/fasttree/FastTree && \
wget http://www.microbesonline.org/fasttree/FastTreeDbl && \
wget http://www.microbesonline.org/fasttree/FastTreeMP && \
chmod +x /usr/local/bin/FastTree*

RUN mkdir -p /working
WORKDIR /working