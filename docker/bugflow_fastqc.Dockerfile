FROM oxfordmmm/bugflow_base:latest  

################## METADATA ###################### 
LABEL about.summary="Image for fastqc prcoesses of BUGflow: nextflow based pipeline for processing bacterial sequencing data"

################## MAINTAINER ###################### 
LABEL maintainer="David Eyre <david.eyre@bdi.ox.ac.uk>"

RUN mamba env update -n base --file bugflow_conda/bugflow_fastqc.yaml
