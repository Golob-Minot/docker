# VERSION               golob/prokka:1.1.14__bcw.0.3.1
FROM      ubuntu:16.04
RUN export DEBIAN_FRONTEND=noninteractive

RUN apt-get update && apt-get install -y \
python \
python-pip \
wget \
libdatetime-perl \
libxml-simple-perl \
libdigest-md5-perl \
bioperl \
pigz \
awscli \
default-jre \
less \
&& apt-get clean \
&& apt-get purge \
&& rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/*

# Install blast+
RUN cd /usr/local/bin && \
wget ftp://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/2.9.0/ncbi-blast-2.9.0+-x64-linux.tar.gz && \
tar xzvf ncbi-blast-2.9.0+-x64-linux.tar.gz && \
ln -s /usr/local/bin/ncbi-blast-2.9.0+/bin/* /usr/local/bin && \
rm ncbi-blast-2.9.0+-x64-linux.tar.gz

# Install Prokka
RUN mkdir /usr/prokka && \
	cd /usr/prokka && \
	wget https://github.com/tseemann/prokka/archive/v1.14.0.tar.gz && \
	tar xzvf v1.14.0.tar.gz && \
    rm -fr prokka-1.14.0/perl5 && mkdir prokka-1.14.0/perl5 && \
	prokka-1.14.0/bin/prokka --setupdb
ENV PATH="/usr/prokka/prokka-1.14.0/bin:${PATH}"

# Install barrnap (for rRNA detection in Prokka)
RUN mkdir -p /src/barrnap/ && cd /src/barrnap/ && \
wget https://github.com/tseemann/barrnap/archive/0.9.tar.gz && \
tar xzvf 0.9.tar.gz && \
cp /src/barrnap/barrnap-0.9/bin/barrnap /usr/local/bin/ && \
mkdir -p /usr/local/db/ && \
cp /src/barrnap/barrnap-0.9//db/* /usr/local/db/ && \
cd / && rm -rf /src/barrnap


# Install tbl2asn
RUN mkdir -p /src/prokka && \
    cd /src/prokka && \
	wget ftp://ftp.ncbi.nlm.nih.gov/toolbox/ncbi_tools/converters/by_program/tbl2asn/linux64.tbl2asn.gz && \
	gunzip linux64.tbl2asn.gz && \
	chmod +x linux64.tbl2asn && \
	mv /src/prokka/linux64.tbl2asn /usr/local/bin/tbl2asn && \
    rm -rf /src/prokka 

# For sciluigi
RUN pip install bucket_command_wrapper==0.3.1

# Various directories
RUN mkdir -p ~/.aws
RUN mkdir -p /mnt/inputs/s3 && mkdir -p /mnt/inputs/file 
RUN mkdir -p /mnt/outputs/s3 && mkdir -p /mnt/outputs/file 
RUN mkdir /scratch && mkdir /working
RUN mkdir /fh/ && mkdir /app