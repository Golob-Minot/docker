# DADA2
#
# VERSION               golob/dada2:1.12.0.ub.1804__bcw.0.3.1

FROM      ubuntu:18.04

env DEBIAN_FRONTEND=noninteractive
RUN apt-get update && apt-get install -y \
python3-dev \
python3-pip \
gnupg \
software-properties-common \
libssl-dev \
libssh2-1-dev \
curl \
libcurl4-gnutls-dev \
libgit2-dev \
libxml2-dev \
libxslt-dev

RUN apt-key adv --keyserver keyserver.ubuntu.com --recv-keys E298A3A825C0D65DFD57CBB651716619E084DAB9 && \
add-apt-repository 'deb https://cloud.r-project.org/bin/linux/ubuntu bionic-cran35/'

RUN apt-get install -y \
r-base=3.6.1-3bionic \
&& apt-get clean \
&& apt-get purge \
&& rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/*


RUN pip3 install \
awscli>=1.15.14 \
boto3>=1.7.14 \
bucket_command_wrapper==0.3.1

# Dada2 install from git
RUN     Rscript -e 'install.packages("devtools")'
RUN     Rscript -e 'install.packages("BiocManager"); BiocManager::install(update=TRUE, ask=FALSE)'
RUN     Rscript -e 'library(devtools); devtools::install_github("benjjneb/dada2", ref="v1.12");'

# For singularity on the hutch cluster
RUN mkdir /fh && mkdir /app && mkdir /src
RUN mkdir -p /mnt/inputs/file && mkdir -p /mnt/outputs/file && mkdir /scratch && mkdir /working