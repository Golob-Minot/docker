# barcodecop
#
# VERSION               0.5.3

FROM      alpine:3.6
RUN mkdir /fh
RUN mkdir /app
RUN apk add --no-cache py2-pip=9.0.1-r1 bash git
RUN pip install pip --upgrade
RUN mkdir /src
WORKDIR /src
RUN git clone https://github.com/jgolob/barcodecop
WORKDIR /src/barcodecop
RUN git checkout csv-out
RUN python setup.py install
RUN apk del git