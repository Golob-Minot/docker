# golob/cutadapt
#
# VERSION               golob/cutadapt:2.3__bcw.0.3.0_al38B_FH_B

FROM biocontainers/biocontainers:v1.1.0_cv2
USER root
# Create some mount points in the container for use by bucket-command-wrapper
RUN mkdir -p /mnt/inputs/file && \
    mkdir -p /mnt/outputs/file && \
    mkdir /scratch && \
    mkdir /working && \
    mkdir -p /mnt/inputs/s3/ && \
    mkdir -p /mnt/outputs/s3/ && \
    mkdir /fh/ && mkdir /app/

# Install at least python3 (used by BCW). It's OK to change the specific version of python3 used. 
RUN apt-get update && \
    apt-get install -y python3 python3-pip python3-dev gcc make && \
    apt-get clean && \
    apt-get purge && \
    rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/*

# Since we are ONLY installing python3 link to it to make it the default python
RUN ln -s /usr/bin/python3 /usr/bin/python

# Install bucket_command_wrapper via pip, along with boto3 / awscli if we want to use AWS at all
RUN pip3 install \
awscli>=1.15.14 \
boto3>=1.7.14 \
bucket_command_wrapper==0.3.0 \
cutadapt==2.3
