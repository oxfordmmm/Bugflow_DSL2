FROM oxfordmmm/bugflow_base:latest  

################## METADATA ###################### 
LABEL about.summary="Image for snippy prcoesses of BUGflow: nextflow based pipeline for processing bacterial sequencing data"

################## MAINTAINER ###################### 
LABEL maintainer="David Eyre <david.eyre@bdi.ox.ac.uk>"

RUN apt-get clean all && \
 apt-get update && \
 apt-get install -y libgsl-dev && \
 apt-get clean && \
 apt-get purge && \
 rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/*

RUN mamba env update -n base --file bugflow_conda/bugflow_snippy.yaml
