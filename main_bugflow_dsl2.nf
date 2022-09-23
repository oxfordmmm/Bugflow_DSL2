#!/usr/bin/env nextflow

/* 
DSL2 version of Bugflow: A pipeline for mapping followed by variant calling and de novo assembly of Illumina short read libraries

QC
 - Fastp
 - FastQC
 - MultiQC
 - Quast

Mapping and Variant calling
 - snippy
 
Assembly
 - shovill (spades) 
*/

// enable DSL2
nextflow.enable.dsl = 2

/*
#==============================================
Modules
#==============================================
*/

include { RAWFASTQC; FASTP; CLEANFASTQC; MULTIQC  } from './modules/processes-bugflow_dsl2.nf'
include { ASSEMBLY } from './modules/processes-bugflow_dsl2.nf'
include { QUAST } from './modules/processes-bugflow_dsl2.nf'
include { SNIPPYFASTQ1 } from './modules/processes-bugflow_dsl2.nf'
include { SNIPPYFASTA } from './modules/processes-bugflow_dsl2.nf' 
include { SNIPPYCORE } from './modules/processes-bugflow_dsl2.nf'

/*
#==============================================
Parameters
#==============================================
*/

params.ref = " "
params.index = " "
params.reads = " "
params.outdir = " "
params.contigs = " "

/*
#==============================================
Channels
#==============================================
*/

Channel.fromPath(params.ref, checkIfExists:true)
       .set{refFasta}
       //.view()

Channel.fromFilePairs(params.reads, checkIfExists: true)
       .map{it}
       .set{reads}

//Channel.fromPath(params.contigs, checkIfExists:true)
       //.set{assembly}
       //.view()

//return

/*
#==============================================
Workflows
#==============================================
*/

workflow shovill {
    
    main:
    RAWFASTQC(reads)
    FASTP(reads)
    CLEANFASTQC(FASTP.out.reads)
    MULTIQC(RAWFASTQC.out.mix(CLEANFASTQC.out).collect())
    ASSEMBLY(FASTP.out.reads)
    QUAST(ASSEMBLY.out)
}

//return

workflow snippy_fastq {
    
    main:
    //FASTP(reads)
    //SNIPPYFASTQ(FASTP.out.reads, refFasta)
    SNIPPYFASTQ1(reads, refFasta)
    //emit:
    //SNIPPYFASTQ1.out // results
}       

workflow snippy_fasta {
    
    main:
    FASTP(reads)
    ASSEMBLY (FASTP.out.reads)
    SNIPPYFASTA(ASSEMBLY.out, refFasta)
    emit:
    SNIPPYFASTA.out // results
}  

//return

workflow  snippy_core {
    FASTP(reads)
    SNIPPYFASTQ(FASTP.out.reads, refFasta)
    SNIPPYCORE(SNIPPYFASTQ.out.collect(), refFasta)        
}

workflow qc_contigs {
    QUAST(assembly)
}
