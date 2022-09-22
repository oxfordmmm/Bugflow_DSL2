#!/usr/bin/env nextflow

// enable DSL2
nextflow.enable.dsl = 2

// initial fastQC
process RAWFASTQC {
	cpus 4

	//conda '/home/ubuntu/miniconda3/envs/fastqc_env'
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


//Read cleaning with Fastp

process FASTP {
	cpus 8

    //conda './fastp_env.yaml'
    //conda '/home/ubuntu/miniconda3/envs/fastp_env'
    

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

// QC clean reads

process CLEANFASTQC {
	cpus 4
	
	//conda '/home/ubuntu/miniconda3/envs/fastqc_env'

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

process MULTIQC {
	
	//conda '/home/ubuntu/miniconda3/envs/multiqc_env'

	tag {"MultiQC raw and clean FastQC files"}

	publishDir "$params.outdir/multiqc/", mode: "symlink"

	input:
    path ('*/*_fastqc/*')
    

    output:
    file "*multiqc_report.html"
    file "*_data"

    script:

    """
    multiqc . -m fastqc
    """
}


process ASSEMBLY {
  	cpus 8

	tag { "assemble ${uuid}" }

	//conda './shovill_env.yaml'
  	//conda '/home/ubuntu/miniconda3/envs/shovill_env'
  
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

process QUAST  {
	tag { " QC assembly using Quast" }
    
    publishDir "$params.outdir/quast", mode: 'symlink'
    
    input:
    path("*_contigs.fa")  
    
    output:
    path("quast")

    script:
    """
    quast --threads ${task.cpus} --output-dir quast *.fa
    """
}

process SNIPPYFASTQ {
	cpus 4

	tag { "call snps from FQs: ${uuid}" }

	//conda '/home/ubuntu/miniconda3/envs/snippy_env'
	    
	publishDir "$params.outdir/snps/", mode: "symlink"

    input:
    tuple val(uuid), path(reads)
    path(refFasta)

    output:
	path("${uuid}_snippy/*") // whole output folder

    """
    snippy --cpus $task.cpus --outdir ${uuid}_snippy --prefix ${uuid} --reference ${refFasta} --R1 ${reads[0]} --R2 ${reads[1]} 
    """

}

process SNIPPYFASTA {
	cpus 4

	tag { "call snps from contigs: ${uuid}" }

	//conda '/home/ubuntu/miniconda3/envs/snippy_env'
        
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

process SNIPPYCORE {
    tag { "create snp alignments: SnippyCore" }

    publishDir "${params.outdir}/snippy_core", mode: "copy", pattern: "snp.core.vcf"
    publishDir "${params.outdir}/snippy_core", mode: "copy", pattern: "snp.core.fasta"
    publishDir "${params.outdir}/snippy_core", mode: "symlink", pattern: "wgs.core.fasta"

    input:
    path(*_snippy)  // collected list of snippy output directories
    path(refFasta)

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



