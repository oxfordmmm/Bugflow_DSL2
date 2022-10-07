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
	
	publishDir "$params.outdir/raw_fastqc", mode: 'copy'

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
    publishDir "$params.outdir/clean_fastqs/", mode: "copy"

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

process MULTIQC_READS {
	
	//conda './conda/multiqc.yaml'
    conda '/home/ubuntu/miniconda3/envs/multiqc_env'

	tag {"Collate and summarize QC files"}

	publishDir "$params.outdir/multiqc_reads/", mode: "copy"

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
  
  	publishDir "$params.outdir/assemblies/", mode: "copy"

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

process QUAST_FROM_READS  {
    
    tag { " QC assembly using Quast" }

    conda './conda/quast.yaml'
    
    publishDir "$params.outdir/quast", mode: 'copy'
    
    input:
    tuple val(uuid), path(assembly)  
    
    output:
    path("quast_${uuid}")

    script:
    """
    quast.py  ${assembly} -o quast_${uuid}
    """
}

process QUAST_FROM_CONTIGS  {
    
    tag { " QC assembly using Quast" }

    conda './conda/quast.yaml'
    
    publishDir "$params.outdir/quast", mode: 'copy'
    
    input:
    path(assembly)  
    
    output:
    path("quast_${assembly}")

    script:
    """
    quast.py --threads ${task.cpus} ${assembly} -o quast_${assembly}
    """
}

process MULTIQC_CONTIGS {
	
	//conda './conda/multiqc.yaml'
    conda '/home/ubuntu/miniconda3/envs/multiqc_env'

	tag {"Collate and summarize QC files"}

	publishDir "$params.outdir/multiqc_contigs/", mode: "copy"

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
#===============================================
AMRG and Plasmid Type Profiling
#===============================================
*/

process AMR_PLM_FROM_READS {
    
    tag { "AMR finding with Abricate" }
    
    conda './conda/abricate.yaml'

    publishDir "$params.outdir/amr_plasmid/", mode: 'copy'
    
    input:
    tuple val(uuid), path(assembly)  
    
    
    output:
    path("${uuid}_card.tab")
    path("${uuid}_resfinder.tab")
    path("${uuid}_plasmidfinder.tab")
    path("${uuid}_card_summary.tsv")
    path("${uuid}_resfinder_summary.tsv")
    path("${uuid}_plasmidfinder_summary.tsv")

    script:
    """
    abricate  --db card ${assembly} > ${uuid}_card.tab
    abricate --summary ${uuid}_card.tab > ${uuid}_card_summary.tsv
    abricate  --db resfinder ${assembly} > ${uuid}_resfinder.tab
    abricate --summary ${uuid}_resfinder.tab > ${uuid}_resfinder_summary.tsv
    abricate  --db plasmidfinder ${assembly} > ${uuid}_plasmidfinder.tab
    abricate --summary ${uuid}_plasmidfinder.tab > ${uuid}_plasmidfinder_summary.tsv
    """
}
process AMR_PLM_FROM_CONTIGS {
    
    tag { "AMR finding with Abricate" }
    
    conda './conda/abricate.yaml'

    publishDir "$params.outdir/amr_plasmid/", mode: 'copy'
    
    input:
    tuple val(uuid), path(assembly)  
    
    
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
Determine MLST
#===============================================
*/

process MLST_FROM_READS {
    cpus 4

    tag {"MLST: ${uuid}"}

    conda 'conda/mlst.yaml'

    publishDir "$params.outdir/mlst/", mode: 'copy'
    
    input:
    tuple val(uuid), path(assembly)  
    
    
    output:
    path("${uuid}_ST.tsv")

    """
    mlst --scheme $params.mlstdb ${assembly} > ${uuid}_ST.tsv
    """

}

process MLST_CDIFF_FROM_READS {
    cpus 4

    tag {"MLST: ${assembly}"}

    conda 'conda/mlst.yaml'

    publishDir "$params.outdir/mlst/", mode: 'copy'
    
    input:
    tuple val(uuid), path(assembly)  
    
    
    output:
    path("${uuid}_ST.tsv")

    """
    mlst --scheme cdifficile ${assembly} > ${uuid}_ST.tsv
    """

}

process MLST_FROM_CONTIGS {
    cpus 4

    tag {"MLST: ${assembly}"}

    conda 'conda/mlst.yaml'

    publishDir "$params.outdir/mlst/", mode: 'copy'
    
    input:
    path(assembly)  
    
    
    output:
    path("${assembly}_ST.tsv")

    """
    mlst --scheme $params.mlstdb ${assembly} > ${assembly}_ST.tsv
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

	//conda './conda/snippy.yaml'
	    
	publishDir "$params.outdir/snps/", mode: "copy"

    input:
    tuple val(uuid), path(reads), path(refFasta)
    

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

	//conda './conda/snippy.yaml'
        
	publishDir "$params.outdir/snps/", mode: "copy"

    input:
    tuple val(uuid), path(fasta)
    path(refFasta)

    output:
    file("${uuid}_snippy/*"), emit: snippy // whole output folder

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

    //conda './conda/snippy.yaml'
    

    publishDir "${params.outdir}/snippy_core", mode: "copy"
    

    input:
    path(snippy)
    path(refFasta)  // collected list of snippy output directories + ref genome
    

    output:
    //path("*")
    path("snp.core.fasta")
    path("snp.core.vcf")
    path("wgs.core.fasta"), emit: for_gubbins

    """
    snippy-core --ref ${refFasta} --prefix core ${params.outdir}/snps/*_snippy/
    mv core.aln snp.core.fasta
    mv core.vcf snp.core.vcf
    snippy-clean_full_aln core.full.aln > wgs.core.fasta

    """

}

process GUBBINS {
    tag { "Remove recombinant segments with Gubbins" }

    
    publishDir "${params.outdir}/gubbins", mode: "copy"
    

    input:
    path(for_gubbins)

    output:
    path("wgs.core.filtered_polymorphic_sites.fasta"), emit: polymorphsites
    path("*")

    script:
    """
    run_gubbins.py ${for_gubbins}

    """
}

process SNP_SITES {

    tag { "Create a SNP-only, non-rec alignment" }

    
    publishDir "${params.outdir}/snp-sites", mode: "copy"
    
    input:
    path(polymorphsites)

    output:
    path("core.full.nonrec.aln"), emit: nonrec

    script:
    """
    snp-sites -c ${polymorphsites} > core.full.nonrec.aln 

    """
    
}

process PHYLOTREE {

    tag { "Create a non-rec tree using FastTree" }

    
    publishDir "${params.outdir}/tree", mode: "copy"
    
    input:
    path(nonrec)

    output:
    path("*")

    script:
    """
    FastTree -nt ${nonrec} > core.full.nonrec.fasttree.tre
    """

}

/*
#==============================================
Index Reference Genome
#==============================================
*/

process INDEXREFERENCE {
    tag {"index reference FASTA"}
    
	input:
    path (refFasta)
	
	output:
	publishDir "$params.outdir"
	path ("${refFasta}_bwa.fai"), emit: bwa_fai
    path("${refFasta}_samtools.fai"), emit: sam_fai
    //path ("*"), emit: ref_index

	script:
	"""
	#bwa index
	bwa index ${refFasta} > ${refFasta}_bwa.fai
	#samtools index
	samtools faidx ${refFasta} > ${refFasta}_samtools.fai
	#blast indexes for self-self blast
	makeblastdb -dbtype nucl -in $refFasta
	#reference mask
    #genRefMask.py -r $refFasta -m 200 -p 95
    #bgzip -c ${refFasta}.rpt.regions > ${refFasta.baseName}.rpt_mask.gz
	#echo '##INFO=<ID=RPT,Number=1,Type=Integer,Description="Flag for variant in repetitive region">' > ${refFasta.baseName}.rpt_mask.hdr
	#tabix -s1 -b2 -e3 ${refFasta.baseName}.rpt_mask.gz
    """
}


/*
#==============================================
Mask Reference Genome
#==============================================
*/



process REFMASK {

	conda '/home/ubuntu/miniconda3/envs/blast_env'
    input:
	path(refFasta)	

	output:
    path("${refFasta.baseName}")
    //path("${refFasta.baseName}.rpt_mask.hdr")
	//path("${refFasta.baseName}.rpt_mask.gz"), emit: refFasta.baseName

	script:
	"""
    /home/ubuntu/Bugflow_DSL2/bin/genRefMask.py -r ${refFasta} -m 200 -p 95
    bgzip -c ${refFasta}.rpt.regions > ${refFasta.baseName}
	echo '##INFO=<ID=RPT,Number=1,Type=Integer,Description="Flag for variant in repetitive region">' > ${refFasta.baseName}.rpt_mask.hdr
	tabix -s1 -b2 -e3 ${refFasta.baseName}.rpt_mask.gz
    """
}

/*
#==============================================
Map reads to Reference genome using BWA
#==============================================
*/

process BWA {
	cpus 8
	
    tag { "map clean ${uuid} reads to reference" }

    publishDir "${params.outdir}/bwa", mode: "copy"
    
	input:
	tuple val(uuid), path(reads), path(bwa_fai)

    output:
    tuple val(uuid), path("${uuid}.aligned.sam"), emit: bwa_mapped
    
	//don't add read group header here results in poorly formatted header
    
	"""
    bwa mem -r 1.5 -O 6 -t ${task.cpus} bwa_fai ${uuid}_clean.1.fq.gz ${uuid}_clean.2.fq.gz | samtools view -uS - | samtools sort > ${uuid}.aligned.sam
    """
}

/*
#==============================================
Remove duplicates using Samtools v.1.9
#==============================================
*/

process REMOVE_DUPLICATES {
    cpus 4

    conda './conda/samtools.yaml'

	tag "remove duplicates ${uuid}"
	
	publishDir "${params.outdir}/bwa", mode: "copy"
    
	input:
    tuple val (uuid), path (bwa_mapped)
    

    output:
    tuple val(uuid), path("${uuid}.bam"), emit: dup_removed

	//sort by name to run fixmate (to remove BWA artefacts) and provide info for markdup
	//sort by position to run markdup (and then remove duplicates)
    """
    samtools sort  -n -o sorted.bam ${uuid}.aligned.sam
    samtools fixmate -m sorted.bam fixed.sorted.bam
    samtools sort -o fixed.resorted.bam fixed.sorted.bam
    samtools markdup -r fixed.resorted.bam ${uuid}.bam
    samtools index ${uuid}.bam
    """
}

/*
#=================================================================
Run Samtools mpileup - creates BCF containing genotype likelihoods 
#=================================================================
*/

process MPILEUP {

    conda './conda/bcftools.yaml'

    publishDir "${params.outdir}/pileup", mode: "copy"

    input:
    tuple val(uuid), path(bam), path(refFasta)
    	
 
    output:
    tuple val(uuid), path("${uuid}.pileup.bcf"), emit: pileup
   
    
    //publishDir "$outputPath/$uuid/bwa_mapped/${refFasta.baseName}/vcf", mode: 'copy'

	//use bcftools mpileup to generate vcf file
	//mpileup genearates the likelihood of each base at each site
 	"""
    bcftools mpileup -Q25 -q30 -E -o40 -e20 -h100 -m2 -F0.002 -Ou -f ${refFasta} ${uuid}.bam > ${uuid}.pileup.bcf
    """

}

/*
#=================================================================
Call SNPs using Samtools call from mpileup file
#=================================================================
*/

process SNP_CALL {
    
    conda './conda/bcftools.yaml'
    

    input:
    tuple val(uuid), path(pileup), path(refFasta)
    
 
    output:
    tuple val(uuid), path("${uuid}.bcf"), path("${uuid}.allsites.bcf"), emit: snps_called
    //path("${uuid}.bcf"), emit: bcf
    //path("${uuid}.allsites.bcf"), emit: allsites
    
    publishDir "${params.outdir}/snps_called_bcf", mode: 'copy'

	//call converts pileup to actual variants in the BCF or VCF file
	//norm, normalises indels
		//-m +any option to join biallelic sites into multiallelic records
		// and do this for any (i.e. SNPs and indels)
	
    script:

    """      
    # call variants only
    # 	-m use multiallelic model
    # 	-v output variants only
    bcftools call --prior 0.01 -Ou -m -v ${uuid}.pileup.bcf | bcftools norm -f ${refFasta} -m +any -Ou -o ${uuid}.bcf
    	
    # call all sites
    bcftools call -Ou -m ${uuid}.pileup.bcf | bcftools norm -f $refFasta -m +any -Ou -o ${uuid}.allsites.bcf
    """

}

/*
#=================================================================
Produce Cleaner SNPs
#=================================================================
*/

process FILTER_SNPS {

    conda './conda/bcftools.yaml'

    publishDir "${params.outdir}/snps_called_vcf", mode: 'copy'

    input:
    //tuple val(uuid), path(bcf),path(allsites) 
    tuple val(uuid), path(snps_called)
    path(refFasta)

 
    output:
    //tuple val(uuid), path("${uuid}.snps.vcf.gz")
    //path("${uuid}.snps.vcf.gz")
    //path("${uuid}.snps.vcf.gz.csi")
    //path("${uuid}.zero_coverage.vcf.gz")
    //path("${uuid}.zero_coverage.vcf.gz.csi"), emit: filtered_snps
    //file ("*")
    path("${uuid}.masked.bcf.gz")
    path("${uuid}.all.vcf.gz")
    path("${uuid}.snps.vcf.gz")
    path("${uuid}.indels.vcf.gz")
    path("${uuid}.zero_coverage.vcf.gz")


	//publishDir "$outputPath/$uuid/bwa_mapped/${refFasta.baseName}/vcf", mode: 'copy', pattern: "${uuid}.snps.*"
	//publishDir "$outputPath/$uuid/bwa_mapped/${refFasta.baseName}/vcf", mode: 'copy', pattern: "${uuid}.indels.*"
	//publishDir "$outputPath/$uuid/bwa_mapped/${refFasta.baseName}/vcf", mode: 'copy', pattern: "${uuid}.zero_coverage.*"
	//publishDir "$outputPath/$uuid/bwa_mapped/${refFasta.baseName}/vcf", mode: 'copy', pattern: "${uuid}.all.*"
	//publishDir "$outputPath/$uuid/bwa_mapped/${refFasta.baseName}/vcf", mode: 'copy'
	
	//use bcftools to filter normalised bcf file of variants from pileup and call
	//use one line for each filter condition and label
	//create index at end for random access and consensus calling
	
	//filters
		// quality >30
		// one read in each direction to support variant
		// not in a repeat region
		// consensus of >90% reads to support alternative allele
		// mask SNPs within 7 bp of INDEL
		// require high quality depth of 5 for call
	
    """
    #annotate vcf file with repetitive regions
	bcftools annotate -a ${refFasta.baseName}.rpt_mask.gz -c CHROM,FROM,TO,RPT -h ${refFasta.baseName}.rpt_mask.hdr ${uuid}.bcf -Ob -o ${uuid}.masked.bcf.gz
    
    #filter vcf
    bcftools filter -s Q30 -e '%QUAL<30' -Ou ${uuid}.masked.bcf.gz | 
        bcftools filter -s HetroZ -e "GT='het'" -m+ -Ou | 
    	bcftools filter -s OneEachWay -e 'DP4[2] == 0 || DP4[3] ==0' -m+ -Ou | 
    	bcftools filter -s RptRegion -e 'RPT=1' -m+ -Ou | 
    	bcftools filter -s Consensus90 -e '((DP4[2]+DP4[3])/(DP4[0]+DP4[1]+DP4[2]+DP4[3]))<=0.9' -m+ -Ou | 
    	bcftools filter -s HQDepth5 -e '(DP4[2]+DP4[3])<=5' -m+ -Oz -o ${uuid}.all.vcf.gz
    
    #create vcf file with just SNPs
    bcftools filter -i 'TYPE="snp"' -m+ -Oz -o ${uuid}.snps.vcf.gz ${uuid}.all.vcf.gz 
    bcftools index ${uuid}.snps.vcf.gz
    
    #create vcf file with just INDELs
    bcftools filter -i 'TYPE!="snp"' -m+ -Oz -o ${uuid}.indels.vcf.gz ${uuid}.all.vcf.gz 
    bcftools index ${uuid}.indels.vcf.gz
    
    #create vcf file with just zero depth sites
    bcftools filter -e 'DP>0' -Oz -o ${uuid}.zero_coverage.vcf.gz ${uuid}.allsites.bcf
    bcftools index ${uuid}.zero_coverage.vcf.gz
    """

}

/*
#=================================================================
Generate a consensus FASTA file
#=================================================================
*/

process CONSENSUS_FA {

        conda './conda/samtools.yaml'
        conda './conda/bcftools.yaml'


	    input:
		tuple val(uuid), path("${uuid}.snps.vcf.gz"), file("${uuid}.snps.vcf.gz.csi"), 
    	file("${uuid}.zero_coverage.vcf.gz"), file("${uuid}.zero_coverage.vcf.gz.csi")
		file(refFasta)
	
	    output:
		path("tmp.bcf.gz")
        path("tmp.fa")
        path("${uuid}.consensus.fa")
	
	    tag "${getShortId(uuid)}"
	    publishDir "$outputPath/$uuid/bwa_mapped/${refFasta.baseName}/fasta", mode: 'copy', pattern: "${uuid}.*"
	    //publishDir "$outputPath/$uuid/bwa_mapped/${refFasta.baseName}/fasta", mode: 'copy'

	    // call consensus sequence
		// -S flag in bcftools filter sets GT (genotype) to missing, with -M flag here
		// setting value to N
	    """
	    #create a temporary bcf file with genotype of filtered variants set to .
	    bcftools filter -S . -e 'FILTER != "PASS"' -Ob -o tmp.bcf.gz ${uuid}.snps.vcf.gz
	    bcftools index tmp.bcf.gz
	
	    #create consensus file with all the sites set to . above replaced as N
	    cat $refFasta | bcftools consensus -H 1 -M "N" tmp.bcf.gz > tmp.fa
	
	    #set all the sites with zero coverage to be -
	    samtools faidx tmp.fa 
	    cat tmp.fa | bcftools consensus -H 1 -M "-" ${uuid}.zero_coverage.vcf.gz > ${uuid}.consensus.fa
	    """
}


process CGMLST_READS_DE {
    tag { "cgMLST Profiling using Hash-cgMLST" }

    // '/home/ubuntu/anaconda3/envs/hash-cgmlst_env'

    publishDir "${params.outdir}/cgmlst", mode: "copy"

    input:
    tuple val(uuid) path(assembly)

    output:
    path("*")

    script:

    """
    /usr/bin/python3 /home/ubuntu/Rev_Bugflow/hash-cgmlst/bin/getCoreGenomeMLST.py -f ${uuid}  -n ${uuid}_hash-cgmlst -s ridom_scheme/files -d ridom_scheme/ridom_scheme.fasta -o ${params.outdir}/cgmlst/${uuid} -b /home/ubuntu/anaconda3/envs/hash-cgmlst_env/bin/blastn
    """
}

process CGMLST_CONTIGS_DE {
    tag { "cgMLST Profiling using Hash-cgMLST" }

    // '/home/ubuntu/anaconda3/envs/hash-cgmlst_env'

    publishDir "${params.outdir}/cgmlst", mode: "copy"

    input:
    path(assembly)

    output:
    path("*")

    script:

    """
    /usr/bin/python3 /home/ubuntu/Rev_Bugflow/hash-cgmlst/bin/getCoreGenomeMLST.py -f ${assembly}  -n ${assembly}_hash-cgmlst -s ridom_scheme/files -d ridom_scheme/ridom_scheme.fasta -o ${params.outdir}/cgmlst/${assembly} -b /home/ubuntu/anaconda3/envs/hash-cgmlst_env/bin/blastn
    """
}

