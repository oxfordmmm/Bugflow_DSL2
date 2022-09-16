#!/usr/bin/env nextflow
// Build indexes for reference fasta file - bwa, samtools, repeatmask


// enable DSL2
nextflow.enable.dsl = 2

//modules
include { RAWFASTQC; FASTP; CLEANFASTQC; MULTIQC  } from './modules/test-bugflow_dsl2.nf'
include { ASSEMBLY } from './modules/test-bugflow_dsl2.nf'
include { SNIPPYFASTQ } from './modules/test-bugflow_dsl2.nf'
include { SNIPPYFASTA } from './modules/test-bugflow_dsl2.nf' 
include { INDEXREFERENCE; BWA} from './modules/test-bugflow_dsl2.nf'
include { SNIPPYCORE } from './modules/snippy_proc.nf'
//parameters
params.ref = "./EC958seq.fasta"
params.index = "work/*/.fai"
params.reads = "/mnt/microbio/HOMES/arund/MERIT/Main/FQs/PE/Clean_PE/*{1,2}.fastq.gz"//"/mnt/microbio/HOMES/arund/MERIT/NF_test/in_fastqs/*{1,2}.fastq.gz" //"./example_data/*{1,2}.fastq.gz"
params.outdir = "./results"
params.fasta = " "

Channel.fromPath(params.ref, checkIfExists:true)
       .set{refFasta}
       //.view()


Channel.fromFilePairs(params.reads, checkIfExists: true)
       .map{it}
       .set{reads}


// workflows
workflow shovill {
    
       main:
       RAWFASTQC(reads)
       FASTP(reads)
       CLEANFASTQC(FASTP.out.reads)
       MULTIQC(RAWFASTQC.out.mix(CLEANFASTQC.out).collect())
       ASSEMBLY(FASTP.out.reads)
}

//return

workflow mapper {
       FASTP(reads)
       INDEXREFERENCE(refFasta)
       BWA(reads, refFasta)

}

//return

workflow snippy_fastq {
    
    main:
        FASTP(reads)
        SNIPPYFASTQ(FASTP.out.reads, refFasta)
    emit:
        SNIPPYFASTQ.out // results
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
        SNIPPYCORE(SNIPPYFASTQ.out.snps.txt.collect(), refFasta)        
}

