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
    label 'fastqc'
    tag {"FastQC raw ${uuid} reads"}
	
	publishDir "$params.outdir/raw_fastqc", mode: 'copy'

	input:
    tuple val(uuid), path(reads) 
	
	output:
	path ("${uuid}_fastqc/${uuid}.txt")
	
	script:
	
	"""
	fastqc --threads ${task.cpus} ${reads[0]}
	fastqc --threads ${task.cpus} ${reads[1]}
	"""

}

process RAWFASTQC_SINGLE {
    label 'fastqc'
	
    tag {"FastQC raw ${uuid} reads"}
	
	publishDir "$params.outdir/raw_fastqc_single", mode: 'copy'

	input:
    tuple val(uuid), path(reads) 
    
	
	output:
    path ("${uuid}_raw_reads_fastqc/${uuid}.txt")
	
	script:
	
	"""
    cat ${reads[0]} ${reads[1]} > ${uuid}_raw_reads.fastq.gz
	fastqc --threads ${task.cpus} ${uuid}_raw_reads.fastq.gz
    unzip ${uuid}_raw_reads_fastqc.zip
    rm ${uuid}_raw_reads_fastqc.zip
    mv ${uuid}_raw_reads_fastqc/fastqc_data.txt ${uuid}_raw_reads_fastqc/${uuid}.txt
	
	"""
}

/*
#==============================================
Read cleaning with Fastp
#==============================================
*/

process FASTP {
    label 'non_perl'
    tag {"filter $uuid reads"}

    input:
    tuple val(uuid), path(reads) 
    
    output:
    tuple val(uuid), path("${uuid}_clean_R*.fastq.gz"), emit: reads
    path("${uuid}.fastp.json"), emit: json
    
    script:
    """
    fastp -i ${reads[0]} -I ${reads[1]} -o ${uuid}_clean_R1.fastq.gz -O ${uuid}_clean_R2.fastq.gz -w ${task.cpus} -j ${uuid}.fastp.json
    """

    stub:
    """
    touch ${uuid}.fastp.json ${uuid}_clean_R1.fastq.gz ${uuid}_clean_R2.fastq.gz fastp.html
    """
}

process FASTP_SINGLE {
    label 'non_perl'

    tag {"filter $uuid reads"}

    input:
    tuple val(uuid), path(reads) 
    
    output:
    tuple val(uuid), path("${uuid}_clean_R*.fastq.gz"), emit: reads
    path("${uuid}.fastp.json"), emit: json
    tuple val(uuid), path("${uuid}.fastq.gz"), emit: cat_fastq

    script:
    """
    fastp -i ${reads[0]} -I ${reads[1]} -o ${uuid}_clean_R1.fastq.gz -O ${uuid}_clean_R2.fastq.gz -w ${task.cpus} -j ${uuid}.fastp.json
    cat ${uuid}_clean_R*.fastq.gz > ${uuid}.fastq.gz
    """

    stub:
    """
    touch ${uuid}.fastp.json  ${uuid}.fastq.gz ${uuid}_clean_R1.fastq.gz ${uuid}_clean_R2.fastq.gz fastp.html
    """
}

//return

/*
#==============================================
QC clean reads
#==============================================
*/

process CLEANFASTQC {
    label 'fastqc'
	tag {"FastQC clean ${uuid} reads"}

    publishDir "$params.outdir/clean_fastqc", mode: 'copy', pattern: "**.txt"

	input:
    tuple val(uuid), path(reads)
	
	output:
	path ("*")
    
	script:
	"""
    fastqc --threads ${task.cpus} ${reads}
	"""

}

process CLEANFASTQC_SINGLE {
    label 'fastqc'
	tag {"FastQC clean ${uuid} reads"}

    publishDir "$params.outdir/clean_fastqc_single", mode: 'copy', pattern: '*.txt'

	input:
    tuple val(uuid), path(cat_fastq)
    
	output:
	path ("*")
    
	script:
	"""
    fastqc --threads ${task.cpus} ${cat_fastq}
    unzip ${uuid}_fastqc.zip
    cp ${uuid}_fastqc/fastqc_data.txt ${uuid}.txt
	"""

    stub:
    """
    touch ${uuid}_fastqc.html ${uuid}_fastqc.zip
    mkdir ${uuid}_fastqc
    touch ${uuid}_fastqc/fastqc.fo ${uuid}_fastqc/fastqc_data.txt ${uuid}_fastqc/fastqc_report.html \
        ${uuid}_fastqc/summary.txt ${uuid}.txt
    """

}

//return

/*
#==============================================
Collate and summarize all read QC files
#==============================================
*/

process MULTIQC_READS {
    label 'non_perl'
	tag {"Collate and summarize QC files"}

	publishDir "$params.outdir/multiqc_reads/", mode: "copy"

	input:
    path ('*')
    
    output:
    path ("*multiqc_report.html")
    path ("*_data")

    script:
    """
    multiqc . 
    """

    stub:
    """
    touch multiqc_report.html
    mkdir multiqc_data
    touch multiqc_data/multiqc.log multiqc_data/multiqc_citations.txt multiqc_data/multiqc_data.json \
        multiqc_data/multiqc_fastqc.txt multiqc_data/multiqc_general_stats.txt multiqc_data/multiqc_sources.txt
    """
}

/*
#==============================================
De novo assembly
#==============================================
*/

process ASSEMBLY {
	tag { "assemble ${uuid}" }
    label 'shovill'

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

    stub:
    """
    touch ${uuid}_contigs.fa
    """
}

/*
#==============================================
QC assembled genomes
#==============================================
*/

process QUAST_FROM_READS {
    
    tag { " QC assembly using Quast" }
    
    label 'quast'
    publishDir "$params.outdir/quast", mode: 'copy', pattern: "*.tsv"
    
    input:
    tuple val(uuid), path(assembly)  
    
    output:
    path("quast_${uuid}"), emit: quast_dir
    path("*.tsv")

    script:
    """
    quast.py  ${assembly} -o quast_${uuid}
    cp quast_${uuid}/report.tsv ${uuid}_Quastreport.tsv
    """

    stub:
    """
    mkdir quast_${uuid}
    touch quast_${uuid}/report.tsv ${uuid}_Quastreport.tsv quast_${uuid}/quast.log quast_${uuid}/report.html \
        quast_${uuid}/report.pdf quast_${uuid}/report.tex quast_${uuid}/report.tsv quast_${uuid}/report.txt \
        quast_${uuid}/transposed_report.tex quast_${uuid}/transposed_report.tsv quast_${uuid}/transposed_report.txt
    """
}

process QUAST_FROM_CONTIGS  {
    
    tag { "QC assembly using Quast" }

    label 'quast'
    publishDir "$params.outdir/quast", mode: 'copy'
    
    input:
    tuple val(uuid), path(assembly)  
    
    output:
    path("quast_${assembly}")

    script:
    """
    quast.py --threads ${task.cpus} ${assembly} -o quast_${assembly}
    """
}

process MULTIQC_CONTIGS {
	
	label 'non_perl'

	tag {"Collate and summarize QC files"}

	input:
    path ('*')
    
    output:
    file "*multiqc_report.html"
    file "*_data"

    script:

    """
    multiqc . 
    """

    stub:
    """
    touch multiqc_report.html
    mkdir multiqc_data
    touch multiqc_data/multiqc.log multiqc_data/multiqc_citations.txt multiqc_data/multiqc_data.json \
        multiqc_data/multiqc_fastqc.txt multiqc_data/multiqc_general_stats.txt multiqc_data/multiqc_sources.txt
    """
}


/*
#===============================================
AMRG and Plasmid Type Profiling
#===============================================
*/

process AMR_PLM_FROM_READS {
    
    tag { "AMR finding with Abricate" }

    label 'abricate'

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
    path("*")

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

    label 'abricate'

    publishDir "$params.outdir/amr_plasmid/", mode: 'copy'
    
    input:
    tuple val(uuid), path(assembly)  
    
    output:
    path("*")

    script:
    """
    abricate  --db card ${assembly} > ${uuid}_card.tab
    abricate --summary ${uuid}_card.tab > ${uuid}_card_amr_summary.tsv
    abricate  --db resfinder ${assembly} > ${uuid}_resfinder_amr.tab
    abricate --summary ${uuid}_resfinder_amr.tab > ${uuid}_resfinder_summary.tsv
    abricate  --db plasmidfinder ${assembly} > ${uuid}_plasmidfinder.tab
    abricate --summary ${uuid}_plasmidfinder.tab > ${uuid}_plasmidfinder_summary.tsv
    """
}

/*
#===============================================
Determine MLST
#===============================================
*/

process MLST_FROM_READS {
    tag {"MLST: ${uuid}"}

    label 'mlst'

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

    tag {"MLST: ${assembly}"}

    label 'mlst'
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

    tag {"MLST: ${assembly}"}

    label 'mlst'
    publishDir "$params.outdir/mlst/", mode: 'copy'
    
    input:
    tuple val(uuid), path(assembly)
    
    
    output:
    path("${uuid}_ST.tsv")

    """
    mlst --scheme $params.mlstdb ${assembly} > ${uuid}_ST.tsv
    """

}

/*
#===============================================
Read mapping and SNP calling from Illumina reads
#===============================================
*/

process SNIPPYFASTQ {
	tag { "call snps from FQs: ${uuid}" }

    label 'snippy'

	publishDir "$params.outdir/snps/", mode: "copy"

    input:
    tuple val(uuid), path(reads), path(refFasta)

    output:
	path("${uuid}_snippy/*"), emit: snippy // whole output folder

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
	tag { "call snps from contigs: ${uuid}" }

    label 'snippy'
        
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
    
    tag { "create snp alignments using snippy-core" }

    label 'snippy'
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
    label 'non_perl'
    tag { "remove recombinant segments with Gubbins" }
    
    publishDir "${params.outdir}/gubbins", mode: "copy"

    input:
    path(for_gubbins)

    output:
    path("wgs.core.filtered_polymorphic_sites.fasta"), emit: polymorphsites
    path("*")

    script:
    """
    run_gubbins.py ${for_gubbins} --threads  $task.cpus --remove-identical-sequences --tree-builder fasttree 
    """
}

/*
#==================================================
Create a SNP-only (non-recombinant) alignment
#==================================================
*/

process SNP_SITES {
    tag { "create a SNP-only, non-rec alignment" }

    label 'non_perl'
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

/*
#==================================================
Create a pairwise SNP matrix
#==================================================
*/

process SNP_DISTS {
    tag { "create a pairwise SNP matrix" }

    label 'non_perl'
    publishDir "${params.outdir}/snp-dists", mode: "copy"

    input:
    path(nonrec)

    output:
    path("*")
    //path("${nonrec}.snp-dists.tsv"), snp_matrix_tsv
    //path("${nonrec}.snp-dists.csv"), snp_matrix_csv
    //path("${nonrec}.snp-dists_molten.csv"), snp_matrix_molten

    script:
    """
    snp-dists ${nonrec} > ${nonrec}.snp-dists.tsv
    snp-dists -c ${nonrec} > ${nonrec}.snp-dists.csv
    snp-dists -c ${nonrec} -m > ${nonrec}.snp-dists_molten.csv
    """
}

/*
#===================================================================
Generate a non-recombinat SNP-based phylogenomic tree using FastTree
#===================================================================
*/

process PHYLOTREE {

    tag { "generate a non-rec tree using FastTree" }

    label 'non_perl'
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
    label 'blast'
    
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
	label 'blast'

    input:
	path(refFasta)	

	output:
    path("${refFasta.baseName}")
    //path("${refFasta.baseName}.rpt_mask.hdr")
	//path("${refFasta.baseName}.rpt_mask.gz"), emit: refFasta.baseName

	script:
	"""
    genRefMask.py -r ${refFasta} -m 200 -p 95
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
    tag { "map clean ${uuid} reads to reference" }
    label 'blast'

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
    tag "remove duplicates ${uuid}"
    label 'blast'
	
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
    label 'blast'

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
    label 'blast'

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
Produce cleaner SNPs
#=================================================================
*/

process FILTER_SNPS {
    label 'blast'

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
    path("${uuid}.masked.bcf.gz"), emit: masked_bcf
    path("${uuid}.all.vcf.gz"), emit: all_vcf
    path("${uuid}.snps.vcf.gz"), emit: snps_vcf
    path("${uuid}.indels.vcf.gz"), emit: indels_vcf
    path("${uuid}.zero_coverage.vcf.gz"), emit: zerocov_vcf


	
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
    label 'blast'

    publishDir "${params.outdir}/consensus_fa", mode: 'copy'

    input:
    tuple val(uuid), path(snps_vcf), path("${uuid}.snps.vcf.gz.csi"), 
    path(zerocov_vcf), file("${uuid}.zero_coverage.vcf.gz.csi")
    file(refFasta)

    output:
    path("tmp.bcf.gz")
    path("tmp.fa")
    path("${uuid}.consensus.fa")

    
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


process HCGMLST_READS_DE {
    tag { "cgMLST Profiling using Hash-cgMLST: ${uuid}" }
    label 'blast'

    publishDir "${params.outdir}/cgmlst", mode: "copy"

    input:
    tuple val(uuid), path(assembly)

    output:
    //path("*")
    path("${uuid}_cgmlst.json"), emit: cgmlst_json
    path("${uuid}_cgmlst.fa"), emit: cgmlst_fasta
    path("${uuid}_cgmlst.profile"), emit: cgmlst_profile

    script:

    """
    getCoreGenomeMLST.py -f ${assembly} -n ${uuid}_hash-cgmlst -s ${params.ridom_scheme_path}/files -d ${params.ridom_scheme_path}/ridom_scheme.fasta -o ${params.outdir}/cgmlst/${uuid} 
    """
}

process HCGMLST_CONTIGS_DE {
    tag { "cgMLST Profiling using Hash-cgMLST" }
    label 'blast'

    publishDir "${params.outdir}/cgmlst", mode: "copy"

    input:
    tuple val(uuid), path(assembly)

    output:
    //path("*")
    path("${uuid}_cgmlst.json"), emit: cgmlst_json
    path("${uuid}_cgmlst.fa"), emit: cgmlst_fasta
    path("${uuid}_cgmlst.profile"), emit: cgmlst_profile

    script:

    """
    getCoreGenomeMLST.py -f ${assembly} -n ${uuid}_hash-cgmlst -s ${params.ridom_scheme_path}/files -d ${params.ridom_scheme_path}/ridom_scheme.fasta -o ${uuid} 
    compareProfiles.py -i ${params.outdir}/cgmlst/ -o hash-cgmlst_profile_comparisons.txt
    compareProfilesExclude.py -i ${params.outdir}/cgmlst/ -o hash-cgmlst_profile_comparisons_exclude26genes.txt
    """
}

process CDIFF_AMRG_BLASTN_READS {
    tag { "annotate ${uuid} contigs with custom Cdiff AMRG db" }
    label 'blast'

    publishDir "${params.outdir}/cdiff_blastn/", mode: "copy"

    input:
    tuple val(uuid), path(reads)
    path(cdiff_amr_fasta)

    output:
    path("cdiffamr-*_blastn.tsv"), emit: blastn

    script:

    """
    makeblastdb -in ${cdiff_amr_fasta} -parse_seqids  -title "C. diff AMRG db" -dbtype nucl -out cdiffamr
    blastn -query ${params.outdir}/assemblies/${uuid}_contigs.fa -db cdiffamr -out cdiffamr-${uuid}.tsv -perc_identity 95 -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore"
    echo -e "qseqid\tsseqid\tpident\tlength\tmismatch\tgapopen\tqstart\tqend\tsstart\tsend\tevalue\tbitscore" > cdiffamr-${uuid}_blastn.tsv && cat cdiffamr-${uuid}.tsv >> cdiffamr-${uuid}_blastn.tsv
    """
    
}


process SUMMARY_BLASTN {
    tag { "cat BLASTN results" }
    label 'blast'

    publishDir "${params.outdir}/cdiff_blastn", mode: "copy"

    input:
    path(blastn)

    output:
    path("*")
    
    script:
    """
    cat ${blastn} > summary_AMRG_blastn_report.tsv
    """
}

process AMR_ABRFORMAT {
    tag { "AMR finding using custom C. diff db" }
    label 'abricate'

    publishDir "$params.outdir/amr_cdiff_abr/", mode: 'copy'
    
    input:
    tuple val(uuid), path(assembly)  
    
    
    output:
    path("*")
    
    script:
    """
    abricate  --db cdiffamr ${assembly} > ${uuid}_cdiffamr.tab
    abricate --summary  *_cdiffamr.tab > cdiffamr_amr_summary.tsv
    """
}

process AMRFINDERPLUS_CDIFF {
    tag { "AMR finding using AMRFinderPlus" }
    label 'amr_finder'

    publishDir "$params.outdir/amr_cdiff_pointmuts/", mode: 'copy'
    
    input:
    tuple val(uuid), path(assembly)  
    
    
    output:
    path("*")

    script:
    """
    amrfinder -u
    amrfinder --organism Clostridioides_difficile -n ${assembly} > ${uuid}_forpointmuts.tsv
    """

}


process PLATON_READS {
    tag { "ID plasmids in short-read draft assemblies" }
    label 'platon'

    publishDir "${params.outdir}/platon", mode: "copy"

    input:
    tuple val(uuid), path(assembly)

    output:
    path("*")

    script:
    """
    platon  ${assembly} --prefix ${uuid} --threads ${task.cpus}
    """ 

}

process PLATON_CONTIGS {
    tag { "ID plasmids in short-read draft assemblies" }
    label 'platon'

    publishDir "${params.outdir}/platon", mode: "copy"

    input:
    tuple val(uuid), path(assembly)

    output:
    path("*")

    script:
    """
    platon  ${assembly} --prefix ${uuid} --threads ${task.cpus}
    """ 

}

process MOBTYPER {
    tag { "MOB-typer and MOB-recon: ${uuid}" }
    label 'mob_suite'

    publishDir "${params.outdir}/mobsuite", mode: "copy"

    input:
    tuple val(uuid), path(assembly)

    output:
    path("*")

    script:
    """
    mob_typer -x -i ${assembly} -o ${uuid}_mobtyper_results.txt
    #mob_recon --infile ${assembly} --outdir ${uuid}_mobrecon_outdir
    """ 

}
