# barcodecop
#
# VERSION               0.4.1A

FROM biocontainers/biocontainers:v1.1.0_cv2
USER root
RUN mkdir /fh && mkdir /app && mkdir /src
RUN mkdir -p /mnt/inputs/file && mkdir -p /mnt/outputs/file && mkdir /scratch && mkdir /working
RUN apt-get update && \
    apt-get install -y python3 python3-pip && \
    apt-get clean && \
    apt-get purge && \
    rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/*
RUN ln -s /usr/bin/python3 /usr/local/bin/python
USER biodocker
RUN pip3 install \
awscli>=1.15.14 \
boto3>=1.7.14 \
bucket_command_wrapper==0.3.0 \
barcodecop==0.4.1
