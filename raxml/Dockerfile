# raxml
#
# VERSION               8.2.11_bcw

FROM      ubuntu:16.04
RUN apt-get update && apt-get install -y wget git build-essential
RUN mkdir /fh && mkdir /app
RUN mkdir /mnt/inputs && mkdir /mnt/outputs
RUN wget https://raw.githubusercontent.com/jgolob/sciluigi/containertask/tools/bucket_command_wrapper.py
RUN chmod +x bucket_command_wrapper.py && cp bucket_command_wrapper.py /usr/bin/
RUN rm bucket_command_wrapper.py
# Get into the actual RAXML stuff
RUN mkdir /raxml && mkdir /raxml/src
WORKDIR /raxml/src
RUN wget https://github.com/stamatak/standard-RAxML/archive/v8.2.11.tar.gz
RUN tar xzvf v8.2.11.tar.gz
WORKDIR /raxml/src/standard-RAxML-8.2.11

# Pthreads versions
RUN make -f Makefile.SSE3.PTHREADS.gcc
RUN cp /raxml/src/standard-RAxML/raxmlHPC-PTHREADS-SSE3 /usr/local/bin
RUN rm *.o && make -f Makefile.AVX.PTHREADS.gcc
RUN cp /raxml/src/standard-RAxML/raxmlHPC-PTHREADS-AVX /usr/local/bin
RUN rm *.o && make -f Makefile.PTHREADS.gcc
RUN cp /raxml/src/standard-RAxML/raxmlHPC-PTHREADS /usr/local/bin
RUN rm *.o && make -f Makefile.AVX.PTHREADS.gcc
RUN cp /raxml/src/standard-RAxML/raxmlHPC-PTHREADS-AVX /usr/local/bin
RUN rm *.o && make -f Makefile.AVX.PTHREADS.gcc
RUN cp /raxml/src/standard-RAxML/raxmlHPC-PTHREADS-AVX /usr/local/bin

# Set the pthreads / avx as the default
RUN ln -s /usr/local/bin/raxmlHPC-PTHREADS-AVX /usr/local/bin/raxml

# cleanup
WORKDIR /
RUN apt-get remove -y wget git build-essential
RUN rm -r /raxml

