FROM oxfordmmm/bugflow_base:latest  

################## METADATA ###################### 
LABEL about.summary="Image for perl elements of BUGflow: nextflow based pipeline for processing bacterial sequencing data"

################## MAINTAINER ###################### 
LABEL maintainer="David Eyre <david.eyre@bdi.ox.ac.uk>"

RUN mamba install -y perl-bioperl>=1.7.2

RUN mamba install -y \
  abricate \
  fastqc \
  mlst \
  quast \
  shovill \
  snippy=4.6 \
  ncbi-amrfinderplus \
  platon \
  mob_suite