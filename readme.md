# Bugflow DSL2
A nextflow pipeline for processing bacterial sequences

## Overview
DSL2 version of Bug-flow DSL1: A pipeline for mapping followed by variant calling, de novo assembly and genome annotation (AMRG, point mutations and cgMLST profile (C. difficile only)) of Illumina short read libraries. The pipeline is developed for use by the Modernising Medical Microbiology consortium based at the University of Oxford.


The workflow has two main parts: de novo based and mapping based. 
- The de novo workflow does read qc (fastp) followed by assembly with Shovill which is qc'd with Quast. MLST and CGMLST are then calculated from the assemblies, and then AMR genes are detected using a mix of AMRFinderPlus and blast (for toxin genes and some specific resistance genes). Kraken2 is also run on the reads to check for contamination.
- The mapping workflow also does read qc (fastp), then maps to the reference using bwa. bcftools is then used to create a pileup, filtering for read and mapping quality. From the pileup bcftools is used to call snps, filter and mask snps, and then create a consensus fasta. 
A check is run for heterozygous calls within the mlst loci as this would indicate a sample with a mix of cdiff strains.
Genome read depth and some quality stats are then calculated. 


The pipeline uses these tools:

QC
 - Fastp v0.23.2
 - FastQC v0.11.9
 - MultiQC v1.12
 - Quast v5.0.2

Mapping and Variant calling
 - BWA mem and bcftools 
 - Snippy v4.6
 
Assembly
 - Shovill v1.1.0 (spades and pilon) 

AMRG Annotation
 - Abricate v0.8
 - AMRFinderPlus v3.11.2
 - BLASTN vs C. difficile curated AMR database

Species Identification
 - Kraken2

## Installation
Requires a local installation of 
* Docker - https://www.docker.com/get-started
* Java version 8 or later (required for nextflow)
* Nextflow - https://www.nextflow.io
* Miniconda - Miniconda3-py37_4.12.0-Linux-x86_64.sh

Tools from Bugflow DSL2 are available from Bioconda.

Install Miniconda locally first,
```
wget https://repo.anaconda.com/miniconda/Miniconda3-py37_4.12.0-Linux-x86_64.sh
bash Miniconda3-py37_4.12.0-Linux-x86_64.sh
```
then create the bug-flow_DSL2 conda environment with all the software needed for running the pipeline.
```
conda create -n Bugflow_DSL2 -c bioconda fastqc fastp multiqc shovill quast mlst abricate snippy samtools bwa bcftools
conda activate Bugflow_DSL2
```

Clone the repository locally
```
git clone https://github.com/aedecano/Bugflow_DSL2
```

Get the docker image
This can be pulled from docker hub
```
docker pull 
```

Alternatively the docker image can be built from the Dockerfile. Within the cloned repository:
```
cd docker
docker build -t .
```
Note the tag has to match in the `nextflow.config` file.

## Running the the pipeline

Bugflow DSL2 has subworkflows for screening and de novo assembly of short reads and accurately calling variants.

To clean and de novo assemble raw Illumina reads:

```
nextflow run main_bugflow_dsl2.nf -entry shovill --reads "[path-to-reads]/*{1,2}.fastq.gz" --outdir "[output_directory]"
```

To call high-quality SNPs from clean reads (customized for C. difficile isolates:):

```
nextflow run main_bugflow_dsl2.nf -entry cdiff_mapping_snpCalling_DE --reads "[path-to-reads]/*{1,2}.fastq.gz" --outdir "[output_directory]" --ref "[you_reference_sequence.fasta]"
```

To call high-quality SNPs from assembled genomes (contigs):

```
nextflow run main_bugflow_dsl2.nf -entry snippy_fasta --reads "[path-to-reads]/*{1,2}.fastq.gz" --outdir "[output_directory]" --ref "[you_reference_sequence.fasta]"
```

## Running the subworkflows on example data

The "example_data" folder included in this repository should contain 2 sets of fastq files. Download these pairs from ENA.

```
cd example_data
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR334/008/SRR3349138/SRR3349138_1.fastq.gz
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR334/008/SRR3349138/SRR3349138_2.fastq.gz
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR334/004/SRR3349174/SRR3349174_1.fastq.gz
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR334/004/SRR3349174/SRR3349174_2.fastq.gz
```
Screen and assemble the example reads

```
nextflow run main.nf -entry cdiff_hcgmlst_amrg_blastn_single --reads "./example_data/*{1,2}.fastq.gz" --outdir "Example_output"
```

Alternatively to resume a partially completed run where the intermediate files have been saved:
```
nextflow run main.nf -entry cdiff_hcgmlst_amrg_blastn_single --reads "./example_data/*{1,2}.fastq.gz" --outdir "Example_output" -resume
```

---
Arun Decano, Jeremy Swann, David Eyre, and Matthew Colpus

arun_decano@ndm.ox.ac.uk
david.eyre@bdi.ox.ac.uk 
crookcs.it@ndm.ox.ac.uk

October 2022
