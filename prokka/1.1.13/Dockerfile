# VERSION               golob/prokka:1.1.13__bcw.0.3.0B
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

# Install the SRA toolkit
RUN cd /usr/local/bin && \
	wget https://ftp-trace.ncbi.nlm.nih.gov/sra/sdk/2.9.6/sratoolkit.2.9.6-ubuntu64.tar.gz && \
	tar xzvf sratoolkit.2.9.6-ubuntu64.tar.gz && \
	ln -s /usr/local/bin/sratoolkit.2.9.6-ubuntu64/bin/* /usr/local/bin/ && \
	rm sratoolkit.2.9.6-ubuntu64.tar.gz

# Install Prokka
RUN mkdir /usr/prokka && \
	cd /usr/prokka && \
	wget https://github.com/tseemann/prokka/archive/v1.13.tar.gz && \
	tar xzvf v1.13.tar.gz && \
    rm -fr prokka-1.13/perl5 && mkdir prokka-1.13/perl5 && \
	prokka-1.13/bin/prokka --setupdb
ENV PATH="/usr/prokka/prokka-1.13/bin:${PATH}"

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
	wget ftp://ftp.ncbi.nih.gov/toolbox/ncbi_tools/converters/by_program/tbl2asn/linux64.tbl2asn.gz && \
	gunzip linux64.tbl2asn.gz && \
	chmod +x linux64.tbl2asn && \
	mv /src/prokka/linux64.tbl2asn /usr/local/bin/tbl2asn && \
    rm -rf /src/prokka 

# For sciluigi
RUN pip install bucket_command_wrapper==0.3.0

# Various directories
RUN mkdir -p ~/.aws
RUN mkdir -p /mnt/inputs/s3 && mkdir -p /mnt/inputs/file 
RUN mkdir -p /mnt/outputs/s3 && mkdir -p /mnt/outputs/file 
RUN mkdir /scratch && mkdir /working
RUN mkdir /fh/ && mkdir /app