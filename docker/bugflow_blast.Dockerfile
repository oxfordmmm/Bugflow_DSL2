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