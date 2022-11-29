FROM oxfordmmm/bugflow_base:latest  

################## METADATA ###################### 
LABEL about.summary="Image for blast processes of BUGflow: nextflow based pipeline for processing bacterial sequencing data"

################## MAINTAINER ###################### 
LABEL maintainer="David Eyre <david.eyre@bdi.ox.ac.uk>"

RUN mamba install -y \
  bwa \
  samtools \
  blast \
  bcftools

RUN git clone https://github.com/davideyre/hash-cgmlst.git && \
    mkdir /bugflow_data && \
    mv /hash-cgmlst/ridom_scheme /bugflow_data && \
    rm -rf /hash-cgmlst