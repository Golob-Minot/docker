# cmalign
#
# VERSION               0.1

FROM      ubuntu:xenial
RUN apt-get update && apt-get install -y wget build-essential
RUN mkdir /fh
RUN mkdir /app
RUN mkdir /cmalign && mkdir /cmalign/src && mkdir /cmalign/data
WORKDIR /cmalign/src
RUN wget http://eddylab.org/infernal/infernal-1.1.2.tar.gz
RUN tar xzvf infernal-1.1.2.tar.gz
WORKDIR /cmalign/src/infernal-1.1.2
RUN ./configure
RUN make && make install
WORKDIR /cmalign/src/infernal-1.1.2/easel
RUN make install
WORKDIR /cmalign/src/infernal-1.1.2/hmmer
RUN make install
RUN wget http://rfam.xfam.org/family/RF00177/cm -O /cmalign/data/SSU_rRNA_bacteria.cm
# Cleanup
WORKDIR /cmalign
RUN rm -R /cmalign/src
RUN apt-get remove -y wget build-essential