
# # bwa
#
# VERSION               golob/bwa:0.7.17__rh.GRCh38__bcw.0.3.0

FROM      golob/bwa:0.7.17__bcw.0.3.0
RUN pip3 install pysam
RUN mkdir -p /db
WORKDIR /db
RUN wget ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/405/GCA_000001405.15_GRCh38/seqs_for_alignment_pipelines.ucsc_ids/GCA_000001405.15_GRCh38_full_plus_hs38d1_analysis_set.fna.bwa_index.tar.gz
RUN tar xzvf GCA_000001405.15_GRCh38_full_plus_hs38d1_analysis_set.fna.bwa_index.tar.gz

