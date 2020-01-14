
RUN cd /src && \
wget https://github.com/fhcrc/deenurp/archive/v0.2.6.tar.gz && \
tar xzvf v0.2.6.tar.gz && \
cd /src/deenurp-0.2.6/ && python setup.py install && \
rm -r /src/*

RUN cd /src && \
wget http://eddylab.org/infernal/infernal-1.1.2-linux-intel-gcc.tar.gz && \
tar xzvf infernal-1.1.2-linux-intel-gcc.tar.gz && \
cp /src/infernal-1.1.2-linux-intel-gcc/binaries/* /usr/local/bin/ && \
rm -r /src/*

RUN cd /usr/local/bin && \
wget http://www.microbesonline.org/fasttree/FastTree && \
wget http://www.microbesonline.org/fasttree/FastTreeDbl && \
wget http://www.microbesonline.org/fasttree/FastTreeMP