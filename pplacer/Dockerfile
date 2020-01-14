# pplacer
#
# VERSION               1.1alpha19

FROM      ubuntu:xenial
RUN mkdir /fh && mkdir /app
RUN apt-get update && apt-get install -y wget unzip python-pip
RUN pip install pip --upgrade
RUN pip install biopython taxtastic
RUN mkdir /pplacer && mkdir /pplacer/src
WORKDIR /pplacer/src
RUN wget https://github.com/matsen/pplacer/releases/download/v1.1.alpha19/pplacer-linux-v1.1.alpha19.zip
RUN unzip pplacer-linux-v1.1.alpha19.zip
RUN cp /pplacer/src/pplacer-Linux-v1.1.alpha19/pplacer /usr/local/bin/  && cp /pplacer/src/pplacer-Linux-v1.1.alpha19/rppr /usr/local/bin/ && cp /pplacer/src/pplacer-Linux-v1.1.alpha19/guppy /usr/local/bin/
WORKDIR /pplacer/src/pplacer-Linux-v1.1.alpha19/scripts
RUN python setup.py install


# cleanup
WORKDIR /
RUN apt-get remove -y wget
RUN rm -r /pplacer
RUN rm -rf /var/lib/apt/lists/*

