# dada2-pplacer
#
# VERSION               0.2

FROM      alpine:3.6
RUN mkdir /fh
RUN mkdir /app
RUN apk add --no-cache --upgrade py2-pip=9.0.1-r1 mercurial=4.3.1-r0 bash gcc python-dev gfortran build-base
RUN pip install pip --upgrade
RUN pip install numpy==1.13.3
RUN pip install pandas==0.20.3
RUN apk del gcc python-dev gfortran build-base
RUN mkdir /scripts
WORKDIR /scripts
RUN hg clone https://bitbucket.org/golob/dada2-pplacer --insecure -u 0.2
RUN cp /scripts/dada2-pplacer/*.py /usr/local/bin/
RUN apk del mercurial
