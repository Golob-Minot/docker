# emirge
#
# VERSION               0.1

FROM      ubuntu
RUN apt-get update && apt-get install -y wget unzip git
RUN mkdir /fh
RUN mkdir /app
RUN mkdir /src
WORKDIR /src
RUN apt-get install -y python-pip
RUN apt-get install -y libncurses5-dev zlib1g-dev libbz2-dev liblzma-dev libcurl4-openssl-dev software-properties-common python-software-properties
RUN apt-add-repository ppa:neufeldlab/ppa && apt-get update && apt-get install -y pandaseq
RUN pip install pip --upgrade
RUN pip install numpy scipy pysam biopython cython
RUN wget https://github.com/torognes/vsearch/releases/download/v2.4.3/vsearch-2.4.3-linux-x86_64.tar.gz
RUN tar xzvf vsearch-2.4.3-linux-x86_64.tar.gz
RUN cp vsearch-2.4.3-linux-x86_64/bin/* /usr/local/bin/
RUN wget https://github.com/samtools/samtools/releases/download/1.5/samtools-1.5.tar.bz2
RUN tar xjvf samtools-1.5.tar.bz2
WORKDIR /src/samtools-1.5
RUN ./configure --prefix=/usr/local/
RUN make; make install
WORKDIR /src/
RUN wget https://sourceforge.net/projects/bowtie-bio/files/bowtie/1.2.1.1/bowtie-1.2.1.1-linux-x86_64.zip/download -O bowtie.zip
RUN unzip bowtie.zip
RUN cp /src/bowtie-1.2.1.1/bowtie* /usr/local/bin/
RUN wget https://golob.org/software/usearch10.0.240_i86linux32
RUN cp /src/usearch10.0.240_i86linux32 /usr/local/bin/usearch
RUN chmod +x /usr/local/bin/usearch
RUN git clone https://github.com/csmiller/EMIRGE
WORKDIR /src/EMIRGE
RUN python setup.py build
RUN python setup.py install


# Cleanup
WORKDIR /
RUN rm -R /src