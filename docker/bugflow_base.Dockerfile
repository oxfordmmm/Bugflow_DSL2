FROM continuumio/miniconda3:latest  

################## METADATA ###################### 
LABEL about.summary="Base image for BUGflow: nextflow based pipeline for processing bacterial sequencing data"

################## MAINTAINER ###################### 
LABEL maintainer="David Eyre <david.eyre@bdi.ox.ac.uk>"

RUN apt-get clean all && \
 apt-get update && \
 apt-get upgrade -y && \
 apt-get install -y file unzip && \
 apt-get clean && \
 apt-get purge && \
 rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/*

RUN conda config --add channels conda-forge && \
    conda config --add channels bioconda
RUN conda install mamba -n base
RUN mamba update conda

COPY bin/* /bugflow_bin/
COPY conda/* /bugflow_conda/
RUN chmod a+x /bugflow_bin/*
ENV PATH="$PATH:/bugflow_bin"