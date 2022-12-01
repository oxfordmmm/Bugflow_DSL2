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

include { RAWFASTQC; RAWFASTQC_SINGLE; FASTP; FASTP_SINGLE; CLEANFASTQC; CLEANFASTQC_SINGLE; MULTIQC_READS; MULTIQC_CONTIGS } from './modules/processes-bugflow_dsl2.nf'
include { ASSEMBLY } from './modules/processes-bugflow_dsl2.nf'
include { QUAST_FROM_READS; QUAST_FROM_CONTIGS } from './modules/processes-bugflow_dsl2.nf'
include { SNIPPYFASTQ } from './modules/processes-bugflow_dsl2.nf'
include { SNIPPYFASTA } from './modules/processes-bugflow_dsl2.nf' 
include { SNIPPYCORE } from './modules/processes-bugflow_dsl2.nf'
include { AMR_PLM_FROM_READS; AMR_PLM_FROM_CONTIGS; PLATON_READS; PLATON_CONTIGS; MOBTYPER; CDIFF_AMRG_BLASTN_READS; SUMMARY_BLASTN; AMR_ABRFORMAT; AMRFINDERPLUS_CDIFF} from './modules/processes-bugflow_dsl2.nf'
include { MLST_FROM_READS; MLST_FROM_CONTIGS; MLST_CDIFF_FROM_READS; HCGMLST_READS_DE; HCGMLST_CONTIGS_DE } from './modules/processes-bugflow_dsl2.nf'
include { INDEXREFERENCE; REFMASK;  BWA; REMOVE_DUPLICATES; MPILEUP; SNP_CALL; FILTER_SNPS; CONSENSUS_FA} from './modules/processes-bugflow_dsl2.nf'
include { GUBBINS; SNP_SITES; SNP_DISTS; PHYLOTREE } from './modules/processes-bugflow_dsl2.nf'

/*
#==============================================
Parameters
#==============================================
*/

params.ref = " "
params.reads = " "
params.outdir = " "
params.contigs = " "
params.mlstdb = "ecoli"
params.prefix = "core"
params.blastn = " "

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
           .map({it -> tuple(it.baseName, it)})
           //.view()
           .set{assembly}
    main:       
    QUAST(assembly)
    MULTIQC(QUAST.out.collect())
}

workflow mlst_amr_plm_reads {
    Channel.fromPath(params.reads, checkIfExists:true)
           //.view()
           .set{reads}
    main:
    FASTP(reads)
    ASSEMBLY(FASTP.out.reads)
    MLST_FROM_READS(ASSEMBLY.out.assembly)
    AMR_PLM_FROM_READS(ASSEMBLY.assembly)
}

workflow mlst_amr_plm_contigs {
    Channel.fromPath(params.contigs, checkIfExists:true)
           .map({it -> tuple(it.baseName, it)})
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
           .map({it -> tuple(it.baseName, it)})
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

workflow  snippy_core_tree {
       Channel.fromFilePairs(params.reads, checkIfExists: true)
           .map{it}
           //.view()
           .set{reads}
       Channel.fromPath(params.ref, checkIfExists:true)
           //.view()       
           .set{refFasta}
       FASTP(reads)
       SNIPPYFASTQ(FASTP.out.reads.combine(refFasta))
       SNIPPYCORE(SNIPPYFASTQ.out, refFasta)
       GUBBINS(SNIPPYCORE.out.for_gubbins)
       SNP_SITES(GUBBINS.out.polymorphsites) 
       SNP_DISTS(SNP_SITES.out.nonrec)
       PHYLOTREE(SNP_SITES.out)      
}



workflow cgmlst_fasta {
       Channel.fromPath(params.contigs, checkIfExists:true)
           .map({it -> tuple(it.baseName, it)})
           //.view()
           .set{assembly}
    main:
    HCGMLST_CONTIGS_DE(assembly)
}

workflow cgmlst_reads {
       Channel.fromFilePairs(params.reads, checkIfExists: true)
           .map{it}
           //.view()
           .set{reads}
       main:
       FASTP(reads)
       ASSEMBLY(FASTP.out.reads)
       HCGMLST_READS_DE(ASSEMBLY.out.assembly)
}

workflow assembly_plmcharac_amr_ns {

       Channel.fromFilePairs(params.reads, checkIfExists: true)
           .map{it}
           //.view()
           .set{reads}
       
       main:
       RAWFASTQC(reads)
       FASTP(reads)
       FASTP_SINGLE(reads)
       //CLEANFASTQC(FASTP.out.reads)
       CLEANFASTQC_SINGLE(FASTP_SINGLE.out.cat_fastq)
       //MULTIQC_READS(RAWFASTQC.out.mix(CLEANFASTQC.out).collect())
       //MULTIQC_READS(CLEANFASTQC.out.collect())
       MULTIQC_READS(CLEANFASTQC_SINGLE.out.collect())
       ASSEMBLY(FASTP.out.reads)
       QUAST_FROM_READS(ASSEMBLY.out.assembly)
       MULTIQC_CONTIGS(QUAST_FROM_READS.out.collect())
       MLST_FROM_READS(ASSEMBLY.out.assembly)
       AMR_PLM_FROM_READS(ASSEMBLY.out.assembly)
       PLATON_READS(ASSEMBLY.out.assembly)
       //MOBTYPER(ASSEMBLY.out.assembly)
       
}


workflow mlst_amr_chrom_plm_contigs_ns {
    Channel.fromPath(params.contigs, checkIfExists:true)
           .map({it -> tuple(it.baseName, it)})
           //.view()
           .set{assembly}

    main:
    MLST_FROM_CONTIGS(assembly)
    AMR_PLM_FROM_CONTIGS(assembly)
    PLATON_CONTIGS(assembly)
    
}
