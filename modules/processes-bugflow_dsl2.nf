#!/usr/bin/env nextflow

/*
#==============================================
Enable DSL2
#==============================================
*/

nextflow.enable.dsl = 2

/*
#==============================================
QC of raw reads
#==============================================
*/

process RAWFASTQC {
	cpus 4

    conda './conda/fastqc.yaml'
	
    tag {"FastQC raw ${uuid} reads"}
	
	publishDir "$params.outdir/raw_fastqc", mode: 'symlink'

	input:
    tuple val(uuid), path(reads) 
	
	output:
	path("*")
	
	script:
	
	"""
	fastqc --threads ${task.cpus} ${reads[0]}
	fastqc --threads ${task.cpus} ${reads[1]}
	"""

}

/*
#==============================================
Read cleaning with Fastp
#==============================================
*/

process FASTP {
	cpus 8

    conda './conda/fastp.yaml'
    

    tag {"filter $uuid reads"}
    publishDir "$params.outdir/clean_fastqs/", mode: "symlink"

    input:
    tuple val(uuid), path(reads) 
    
    output:
    tuple val(uuid), path("${uuid}_clean_R*.fastq.gz"), emit: reads
    path("${uuid}.fastp.json"), emit: json
 
    script:
    """
    fastp -i ${reads[0]} -I ${reads[1]} -o ${uuid}_clean_R1.fastq.gz -O ${uuid}_clean_R2.fastq.gz -w ${task.cpus} -j ${uuid}.fastp.json
    
    """
}

//return

/*
#==============================================
QC clean reads
#==============================================
*/

process CLEANFASTQC {
	cpus 4
	
    conda './conda/fastqc.yaml'
	
	tag {"FastQC clean ${uuid} reads"}

	input:
    tuple val(uuid), path(reads)
	
	output:
	path ("*")
	
	publishDir "$params.outdir/clean_fastqc", mode: 'copy', pattern: "${uuid}*"
	
	"""
	fastqc --threads ${task.cpus} ${reads}
	"""

}

//return

/*
#==============================================
Collate and summarize all read QC files
#==============================================
*/

process MULTIQC {
	
	conda './bin/multiqc.yaml'

	tag {"Collate and summarize QC files"}

	publishDir "$params.outdir/multiqc/", mode: "symlink"

	input:
    path ('*')
    
    output:
    file "*multiqc_report.html"
    file "*_data"

    script:

    """
    multiqc . 
    """
}

/*
#==============================================
De novo assembly
#==============================================
*/

process ASSEMBLY {
  	cpus 8

	tag { "assemble ${uuid}" }

	conda './conda/shovill.yaml'
  
  	publishDir "$params.outdir/assemblies/", mode: "symlink"

  	input:
  	tuple val(uuid), path(reads)

  	output:
  	tuple val(uuid), path("${uuid}_contigs.fa"), emit: assembly

  	script:
 	"""
  	shovill --R1 ${reads[0]} --R2 ${reads[1]} --outdir shovill_${uuid} --cpus ${task.cpus}
	mv shovill_${uuid}/contigs.fa ${uuid}_contigs.fa
  	"""
}

/*
#==============================================
QC assembled genomes
#==============================================
*/

process QUAST  {
    
    tag { " QC assembly using Quast" }

    conda './conda/quast.yaml'
    
    publishDir "$params.outdir/quast", mode: 'symlink'
    
    input:
    path(assembly)  
    
    output:
    path("quast_${assembly}")

    script:
    """
    quast.py --threads ${task.cpus} ${assembly} -o quast_${assembly}
    """
}

/*
#===============================================
AMRG and Plasmid Type Profiling
#===============================================
*/

process AMR_PLM {
    
    tag { "AMR finding with Abricate" }
    
    conda './conda/abricate.yaml'

    publishDir "$params.outdir/abricate/", mode: 'symlink'
    
    input:
    path(assembly)  
    
    
    output:
    path("${assembly}_card.tab")
    path("${assembly}_resfinder.tab")
    path("${assembly}_plasmidfinder.tab")
    path("${assembly}_card_summary.tsv")
    path("${assembly}_resfinder_summary.tsv")
    path("${assembly}_plasmidfinder_summary.tsv")

    script:
    """
    abricate  --db card ${assembly} > ${assembly}_card.tab
    abricate --summary ${assembly}_card.tab > ${assembly}_card_summary.tsv
    abricate  --db resfinder ${assembly} > ${assembly}_resfinder.tab
    abricate --summary ${assembly}_resfinder.tab > ${assembly}_resfinder_summary.tsv
    abricate  --db plasmidfinder ${assembly} > ${assembly}_plasmidfinder.tab
    abricate --summary ${assembly}_plasmidfinder.tab > ${assembly}_plasmidfinder_summary.tsv
    """
}




/*
#===============================================
Read mapping and SNP calling from Illumina reads
#===============================================
*/

process SNIPPYFASTQ {
	cpus 4

	tag { "call snps from FQs: ${uuid}" }

	conda './conda/snippy.yaml'
	    
	publishDir "$params.outdir/snps/", mode: "symlink"

    input:
    tuple val(uuid), path(reads), path(refFasta)
    //path(refFasta)

    output:
	path("${uuid}_snippy/*") // whole output folder

    """
    snippy --cpus $task.cpus --outdir ${uuid}_snippy --prefix ${uuid} --reference ${refFasta} --R1 ${reads[0]} --R2 ${reads[1]} 
    """

}


/*
#==================================================
Read mapping and SNP calling from assembled contigs
#==================================================
*/

process SNIPPYFASTA {
	cpus 4

	tag { "call snps from contigs: ${uuid}" }

	conda './conda/snippy.yaml'
        
	publishDir "$params.outdir/snps/", mode: "symlink"

    input:
    tuple val(uuid), path(fasta)
    path(refFasta)

    output:
    file("${uuid}_snippy/*") // whole output folder

    """
    snippy --cpus $task.cpus --report --outdir $uuid --prefix $uuid --reference ${refFasta} --ctgs $fasta 
    """
}

/*
#==================================================
SNP alignment
#==================================================
*/

process SNIPPYCORE {
    tag { "create snp alignments: SnippyCore" }

    conda './conda/snippy.yaml'

    publishDir "${params.outdir}/snippy_core", mode: "copy", pattern: "snp.core.vcf"
    publishDir "${params.outdir}/snippy_core", mode: "copy", pattern: "snp.core.fasta"
    publishDir "${params.outdir}/snippy_core", mode: "symlink", pattern: "wgs.core.fasta"

    input:
    path("*"), path(refFasta)  // collected list of snippy output directories + ref genome
    

    output:
    path("wgs.core.fasta")
    path("snp.core.fasta")
    path("snp.core.vcf")

    """
    snippy-core --ref ${refFasta} --prefix core *_snippy
    mv core.aln snp.core.fasta
    mv core.vcf snp.core.vcf
    snippy-clean_full_aln core.full.aln > wgs.core.fasta

    """

}



