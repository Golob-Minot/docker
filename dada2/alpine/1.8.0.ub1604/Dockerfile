# DADA2
#
# VERSION               1.8.0__bcw.0.3.0U

FROM      ubuntu:18.04
# For singularity on the hutch cluster
RUN mkdir /fh && mkdir /app && mkdir /src
RUN mkdir -p /mnt/inputs/file && mkdir -p /mnt/outputs/file && mkdir /scratch && mkdir /working
RUN apt-get update && apt-get install python3-dev R-dev
RUN apk add --no-cache bash \
python3==3.6.6-r0 \
gfortran==6.4.0-r8 \
readline-dev==7.0.003-r0 \
zlib-dev==1.2.11-r1 \
bzip2-dev==1.0.6-r6 \
xz-dev \
pcre-dev \
libcurl \
curl-dev \
build-base
RUN pip3 install \
awscli>=1.15.14 \
boto3>=1.7.14 \
bucket_command_wrapper==0.3.0 

# Install R from source
RUN     mkdir /src/
WORKDIR /src/
RUN     wget http://cran.mtu.edu/src/base/R-3/R-3.5.1.tar.gz
RUN     tar xzvf R-3.5.1.tar.gz
WORKDIR /src/R-3.5.1
RUN ./configure --with-x=no
RUN make
RUN make install


RUN     Rscript -e 'install.packages("devtools", repos="https://cran.mtu.edu")'
RUN     Rscript -e 'library(devtools); devtools::install_github("benjjneb/dada2", ref="v1.8");'

# Manually install dependencies :<
# 'Rcpp' ‘Biostrings’, ‘ggplot2’, ‘data.table’, ‘reshape2’, ‘ShortRead’, ‘RcppParallel’
RUN     Rscript -e "install.packages('Rcpp', repos='https://cran.mtu.edu/')"
RUN     Rscript -e "install.packages('ggplot2', repos='https://cran.mtu.edu/')"
RUN     Rscript -e "source('https://bioconductor.org/biocLite.R'); biocLite('Biostrings')"
RUN     Rscript -e "install.packages('data.table', repos='https://cran.mtu.edu/')"
RUN     Rscript -e "install.packages('reshape2', repos='https://cran.mtu.edu/')"
RUN     Rscript -e "source('https://bioconductor.org/biocLite.R'); biocLite('ShortRead')"
RUN     Rscript -e "install.packages('RcppParallel', repos='https://cran.mtu.edu/')"

# Download dada2 tgz directly so we can hope to version 
RUN wget https://bioconductor.org/packages/3.7/bioc/src/contrib/dada2_1.8.0.tar.gz
RUN R CMD INSTALL dada2_1.8.0.tar.gz

