#!/usr/bin/env nextflow

/* 
DSL2 version of Bugflow: A pipeline for mapping followed by variant calling and de novo assembly of Illumina short read libraries

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
 - BLASTN vs C. difficile curated AMR database (also contains C. diff toxin genes)
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
include { AMR_PLM_FROM_READS; AMR_PLM_FROM_CONTIGS; PLATON_READS; PLATON_CONTIGS; MOBTYPER; CDIFF_AMRG_BLASTN_READS; SUMMARY_BLASTN; AMR_ABRFORMAT; AMRFINDERPLUS_CDIFF} from './modules/processes-bugflow_dsl2.nf'
include { MLST_FROM_READS; MLST_FROM_CONTIGS; MLST_CDIFF_FROM_READS; HCGMLST_READS_DE; HCGMLST_CONTIGS_DE } from './modules/processes-bugflow_dsl2.nf'
include { INDEXREFERENCE; REFMASK;  BWA; REMOVE_DUPLICATES; MPILEUP; SNP_CALL; FILTER_SNPS; CONSENSUS_FA} from './modules/processes-bugflow_dsl2.nf'
include { SNIPPYFASTQ; SNIPPYCORE } from './modules/processes-bugflow_dsl2.nf'
/*
#==============================================
Parameters
#==============================================
*/

params.ref = "/mnt/arun_in_bucket/REFSEQ_Cdiff/Cdiff_630_GCA_000009205.1.fasta"
params.reads = "/mnt/arun_in_bucket/UKHSA_test/*{1,2}_001.fastq.gz"
params.outdir = " "
params.contigs = " "
params.mlstdb = "cdifficile"
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

workflow cdiff_mapping_snpCalling_DE {
       Channel.fromFilePairs(params.reads, checkIfExists: true)
           .map{it}
           //.view()
           .set{reads}
       Channel.fromPath(params.ref, checkIfExists:true)
           //.view()       
           .first()
           .set{refFasta}
       main:
       INDEXREFERENCE(refFasta)     
       REFMASK(refFasta)
       FASTP(reads)
       BWA(FASTP.out.reads.combine(refFasta))
       REMOVE_DUPLICATES(BWA.out)
       MPILEUP(REMOVE_DUPLICATES.out.dup_removed.combine(refFasta))
       SNP_CALL(MPILEUP.out.pileup.combine(refFasta))
       FILTER_SNPS(SNP_CALL.out.snps_called, refFasta, REFMASK.out.masked_ref, REFMASK.out.masked_ref_hdr)
       //FILTER_SNPS.out.filtered_snps.view()
       CONSENSUS_FA(FILTER_SNPS.out.filtered_snps, refFasta)
       //MSA(CONSENSUS_FA.out.collect())
       
}



workflow cdiff_asssembly_mlst_amr_plm {
       Channel.fromFilePairs(params.reads, checkIfExists: true)
           .map{it}
           //.view()
           .set{reads}
       //Channel.fromPath(params.ref, checkIfExists:true)
           //.view()       
           //.set{refFasta}
       main:
       RAWFASTQC(reads)
       FASTP(reads)
       CLEANFASTQC(FASTP.out.reads)
       MULTIQC_READS(RAWFASTQC.out.mix(CLEANFASTQC.out).collect())
       ASSEMBLY(FASTP.out.reads)
       //QUAST_FROM_READS(ASSEMBLY.out.assembly)
       //MULTIQC_CONTIGS(QUAST_FROM_READS.out.collect())
       //CGMLST_READS_DE(ASSEMBLY.out.assembly)
       MLST_CDIFF_FROM_READS(ASSEMBLY.out.assembly)
       AMR_PLM_FROM_READS(ASSEMBLY.out.assembly)
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

workflow cdiff_hcgmlst_amrg_blastn {
       Channel.fromFilePairs(params.reads, checkIfExists: true)
           .map{it}
           //.view()
           .set{reads}
       
       main:
       RAWFASTQC(reads)
       FASTP(reads)
       CLEANFASTQC(FASTP.out.reads)
       //MULTIQC_READS(RAWFASTQC.out.mix(CLEANFASTQC.out).collect())
       MULTIQC_READS(CLEANFASTQC.out.collect())
       ASSEMBLY(FASTP.out.reads)
       HCGMLST_CONTIGS_DE(ASSEMBLY.out.assembly)
       CDIFF_AMRG_BLASTN_READS(ASSEMBLY.out.assembly, params.cdiff_amr_fasta)
}

workflow cdiff_hcgmlst_amrg_blastn_single {
       Channel.fromFilePairs(params.reads, checkIfExists: true)
           .map{it}
           //.view()
           .set{reads}
       
       main:
       RAWFASTQC_SINGLE(reads)
       FASTP(reads)
       FASTP_SINGLE(reads)
       CLEANFASTQC_SINGLE(FASTP_SINGLE.out.cat_fastq)
       //MULTIQC_READS(RAWFASTQC.out.mix(CLEANFASTQC.out).collect())
       //MULTIQC_READS(CLEANFASTQC.out.collect())
       MULTIQC_READS(CLEANFASTQC_SINGLE.out.collect())
       ASSEMBLY(FASTP.out.reads)
       QUAST_FROM_READS(ASSEMBLY.out.assembly)
       MULTIQC_CONTIGS(QUAST_FROM_READS.out.quast_dir.collect())
       MLST_CDIFF_FROM_READS(ASSEMBLY.out.assembly)
       HCGMLST_CONTIGS_DE(ASSEMBLY.out.assembly)
       CDIFF_AMRG_BLASTN_READS(ASSEMBLY.out.assembly, params.cdiff_amr_fasta)
       AMRFINDERPLUS_CDIFF(ASSEMBLY.out.assembly)
       //AMR_ABRFORMAT(ASSEMBLY.out.assembly)
}
