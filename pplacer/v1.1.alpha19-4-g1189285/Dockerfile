# pplacer
#
# VERSION               1.1alpha19rc_BCW_0.3.0B recompiled

FROM      alpine:3.7
RUN mkdir /fh && mkdir /app && mkdir /src
RUN mkdir -p /mnt/inputs/file && mkdir -p /mnt/outputs/file && mkdir /scratch
RUN apk add --no-cache  bash \
python3==3.6.5-r0 \
python3-dev==3.6.5-r0 \
python2 \
unzip==6.0-r3 \
libpq \
build-base==0.5-r0 \
postgresql-bdr-dev==9.4.10_p1-r3 \
perl-dev \
perl \
pigz \
sqlite>=3.21.0-r1

RUN ln -s /usr/bin/python3 /usr/local/bin/python
RUN pip3 install pip --upgrade && pip install wheel
RUN pip3 install \
numpy==1.14.3 \
biopython>=1.68

RUN pip3 install \
psycopg2>=2.7.4 \
awscli>=1.15.14 \
boto3>=1.7.14 \
bucket_command_wrapper==0.3.0 \
taxtastic>=0.8.5 \
csvkit==1.0.3

COPY pplacer-linux-v1.1.alpha19-4-g1189285.zip /src/
WORKDIR /src
RUN unzip pplacer-linux-v1.1.alpha19-4-g1189285.zip
RUN cp guppy /usr/local/bin && cp pplacer /usr/local/bin && cp rppr /usr/local/bin
WORKDIR /src/scripts
RUN python2 setup.py install

RUN mkdir -p /pplacer/refpkg && chmod a+w /pplacer/refpkg
WORKDIR /scratch
RUN rm -r /src/
