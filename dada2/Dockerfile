# DADA2
#
# VERSION               1.4.0

FROM      ubuntu:xenial
# For singularity on the hutch cluster
RUN mkdir /fh
RUN mkdir /app
# Add updated R repositories from CRAN (and appropriate validation key)
RUN       echo "deb http://cran.cnr.berkeley.edu/bin/linux/ubuntu xenial/" >> /etc/apt/sources.list && gpg --keyserver keyserver.ubuntu.com --recv-key E084DAB9 && gpg -a --export E084DAB9 |  apt-key add -
RUN       apt-get update # Update
RUN       apt-get install -y r-base-dev=3.4.2-1xenial1 libcurl4-openssl-dev # install r-base-dev
# Install DADA2 via biocLite (would like to do versioning)
RUN     Rscript -e "source('https://bioconductor.org/biocLite.R'); biocLite('dada2')"
