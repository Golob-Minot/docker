# barcodecop
#
# VERSION               0.4.1

FROM      alpine:3.7
RUN mkdir /fh && mkdir /app && mkdir /src
RUN mkdir -p /mnt/inputs/file && mkdir -p /mnt/outputs/file && mkdir /scratch && mkdir /working
RUN apk add --no-cache  bash \
python3==3.6.3-r9
RUN ln -s /usr/bin/python3 /usr/local/bin/python
RUN pip3 install \
awscli>=1.15.14 \
boto3>=1.7.14 \
bucket_command_wrapper==0.3.0 \
barcodecop==0.4.1
