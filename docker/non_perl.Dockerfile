FROM continuumio/miniconda3:latest  

RUN apt-get clean all && \
 apt-get update && \
 apt-get upgrade -y && \
 apt-get clean && \
 apt-get purge && \
 rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/*

RUN conda config --add channels conda-forge && \
    conda config --add channels bioconda
RUN conda install mamba -n base
RUN mamba update conda
RUN mamba install -y \
  fastp \
  fasttree \
  multiqc \
  snp-dists \
  snp-sites