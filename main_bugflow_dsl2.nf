#!/usr/bin/env nextflow

/* 
DSL2 version of Bugflow: A pipeline for mapping followed by variant calling and de novo assembly of Illumina short read libraries

QC
 - Fastp
 - FastQC
 - MultiQC
 - Quast

Mapping and Variant calling
 - Snippy
 
Assembly
 - Shovill (spades) 
*/

/*
#==============================================
Enable DSL2
#==============================================
*/

nextflow.enable.dsl = 2

/*
#==============================================
Modules
#==============================================
*/

include { RAWFASTQC; FASTP; CLEANFASTQC; MULTIQC_READS; MULTIQC_CONTIGS } from './modules/processes-bugflow_dsl2.nf'
include { ASSEMBLY } from './modules/processes-bugflow_dsl2.nf'
include { QUAST_FROM_READS; QUAST_FROM_CONTIGS } from './modules/processes-bugflow_dsl2.nf'
include { SNIPPYFASTQ } from './modules/processes-bugflow_dsl2.nf'
include { SNIPPYFASTA } from './modules/processes-bugflow_dsl2.nf' 
include { SNIPPYCORE } from './modules/processes-bugflow_dsl2.nf'
include { AMR_PLM_FROM_READS; AMR_PLM_FROM_CONTIGS } from './modules/processes-bugflow_dsl2.nf'
include { MLST_FROM_READS; MLST_FROM_CONTIGS; MLST_CDIFF_FROM_READS; CGMLST_CONTIGS_DE } from './modules/processes-bugflow_dsl2.nf'
include { INDEXREFERENCE; REFMASK;  BWA; REMOVE_DUPLICATES; MPILEUP; SNP_CALL; FILTER_SNPS; CONSENSUS_FA} from './modules/processes-bugflow_dsl2.nf'

/*
#==============================================
Parameters
#==============================================
*/

params.ref = " "
params.reads = " "
params.outdir = " "
params.contigs = " "
params.mlstdb = " "
params.tracedir = " "

/*
#==============================================
Channels
#==============================================
*/

//Channel.fromFilePairs(params.reads, checkIfExists: true)
       //.map{it}
       //.view()
       //.set{reads}

//Channel.fromPath(params.ref, checkIfExists:true)
        //.view()       
       //.set{refFasta}

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
    Channel.fromFilePairs(params.reads, checkIfExists: true)
           .map{it}
           //.view()
           .set{reads}
    
    main:
    RAWFASTQC(reads)
    FASTP(reads)
    CLEANFASTQC(FASTP.out.reads)
    MULTIQC_READS(RAWFASTQC.out.mix(CLEANFASTQC.out).collect())
    ASSEMBLY(FASTP.out.reads)
    QUAST_FROM_READS(ASSEMBLY.out.collect())
    MULTIQC_CONTIGS(QUAST_FROM_READS.out.collect())
}

//return

workflow qc_contigs {
    Channel.fromPath(params.contigs, checkIfExists:true)
           //.view()
           .set{assembly}
    main:       
    QUAST(assembly)
    MULTIQC(QUAST.out.collect())
}

workflow mlst_amr_plm_reads {
    Channel.fromPath(params.contigs, checkIfExists:true)
           //.view()
           .set{assembly}
    main:
    FASTP(reads)
    ASSEMBLY(FASTP.out.reads)
    MLST_FROM_READS(ASSEMBLY.out.assembly)
    AMR_PLM_FROM_READS(ASSEMBLY.assembly)
}

workflow mlst_amr_plm_contigs {
    Channel.fromPath(params.contigs, checkIfExists:true)
           //.view()
           .set{assembly}
    main:
    MLST_FROM_CONTIGS(assembly)
    AMR_PLM_FROM_CONTIGS(assembly)
    
}

workflow snippy_fastq {
    Channel.fromFilePairs(params.reads, checkIfExists: true)
           .map{it}
           //.view()
           .set{reads}
    
    Channel.fromPath(params.ref, checkIfExists:true)
           //.view()       
           .set{refFasta}
    
    main:
    FASTP(reads)
    SNIPPYFASTQ(FASTP.out.reads.combine(refFasta))
    SNIPPYCORE(SNIPPYFASTQ.out,refFasta) 
}       

workflow snippy_fasta {
    Channel.fromPath(params.contigs, checkIfExists:true)
           //.view()
           .set{assembly}

    Channel.fromPath(params.ref, checkIfExists:true)
           //.view()       
           .set{refFasta}

    main:
    SNIPPYFASTA(assembly.combine(refFasta))
    
    emit:
    SNIPPYFASTA.out // results
}  


//return

workflow  snippy_core {
    FASTP(reads)
    SNIPPYFASTQ(FASTP.out.reads, refFasta)
    SNIPPYCORE(SNIPPYFASTQ.out.combine(refFasta))        
}

workflow cdiff_mapping_snpCalling_DE {
       //Channel.fromFilePairs(params.reads, checkIfExists: true)
           //.map{it}
           //.view()
           //.set{reads}
       Channel.fromPath(params.ref, checkIfExists:true)
           //.view()       
           .set{refFasta}
       main:
       INDEXREFERENCE(refFasta)     
       //REFMASK(refFasta)
       //FASTP(reads)
       //BWA(FASTP.out.reads, REFMASK.out)
       //REMOVE_DUPLICATES(BWA.out)
       //MPILEUP(REMOVE_DUPLICATES.out.dup_removed, refFasta)
       //SNP_CALL(MPILEUP.out.pileup, refFasta)
       //FILTER_SNPS(SNP_CALL.out.snps_called)
       //CONSENSUS_FA(FILTER_SNPS.out.filtered_snps, refFasta)
     
}


workflow cdiff_asssembly_mlst_amr_plm {
       Channel.fromFilePairs(params.reads, checkIfExists: true)
           .map{it}
           //.view()
           .set{reads}
       Channel.fromPath(params.ref, checkIfExists:true)
           //.view()       
           .set{refFasta}
       main:
       RAWFASTQC(reads)
       FASTP(reads)
       CLEANFASTQC(FASTP.out.reads)
       MULTIQC_READS(RAWFASTQC.out.mix(CLEANFASTQC.out).collect())
       ASSEMBLY(FASTP.out.reads)
       //QUAST_FROM_READS(ASSEMBLY.out.assembly)
       //MULTIQC_CONTIGS(QUAST_FROM_READS.out.collect())
       MLST_CDIFF_FROM_READS(ASSEMBLY.out.assembly)
       AMR_PLM_FROM_READS(ASSEMBLY.out.assembly)
}

workflow cgmlst_fasta {
    Channel.fromPath(params.contigs, checkIfExists:true)
           //.view()
           .set{assembly}
    main:
    CGMLST_CONTIGS_DE(assembly)
}
