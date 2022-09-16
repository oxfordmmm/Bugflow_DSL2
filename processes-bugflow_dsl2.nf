#!/usr/bin/env nextflow

// enable DSL2
nextflow.enable.dsl = 2


// initial fastQC
process RAWFASTQC {
	conda '/home/ubuntu/miniconda3/envs/fastqc_env'
	cpus 4
	tag {"FastQC raw ${uuid} reads"}
	
	publishDir "$params.outdir/raw_fastqc", mode: 'symlink'

	input:
        tuple val(uuid), path(reads) 
	
	output:
		path("*")
	
	script:
	
	"""
	cat ${reads[0]} ${reads[1]} > ${uuid}_raw.fq.gz
	fastqc --threads ${task.cpus} ${reads} 
	rm ${uuid}_raw.fq.gz
	"""

}

//return

//Read cleaning with Fastp

process FASTP {
    //conda './fastp_env.yaml'
    conda '/home/ubuntu/miniconda3/envs/fastp_env'
    cpus 8

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

// repeat fastQC
process CLEANFASTQC {
	cpus 4
	conda '/home/ubuntu/miniconda3/envs/fastqc_env'

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
	conda '/home/ubuntu/miniconda3/envs/multiqc_env'

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
  //conda './shovill_env.yaml'
  conda '/home/ubuntu/miniconda3/envs/shovill_env'
  
  cpus 8

  tag { "assemble ${uuid}" }

  publishDir "$params.outdir/assemblies/", mode: "symlink"

  input:
  //path(reads)
  tuple val(uuid), path(reads)

  output:
  tuple val(uuid), path("${uuid}_contigs.fa"), emit: assembly

  script:
  """
  shovill --R1 ${reads[0]} --R2 ${reads[1]} --outdir shovill_${uuid} --cpus ${task.cpus}
  
  mv shovill_${uuid}/contigs.fa ${uuid}_contigs.fa
  """
}

process SNIPPYFASTQ {
		conda '/home/ubuntu/miniconda3/envs/snippy_env'

		cpus 4
		tag { "call snps from FQs: ${uuid}" }
                
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

		conda '/home/ubuntu/miniconda3/envs/snippy_env'
    
		cpus 4
        
        tag { "call snps from contigs: ${uuid}" }

		publishDir "$params.outdir/snps/", mode: "symlink"

        input:
        tuple val(uuid), path(fasta)
        path(refFasta)

        output:
        file("$id") // whole output folder

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


// Map reads to reference genome with BWA MEM

//return

process INDEXREFERENCE {
    tag {"index reference FASTA"}
    
	input:
  	path (refFasta)
	
	output:
	publishDir "$params.outdir"
	path ("*"), emit: ref_index

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
  // reference mask
process refMask {
	input:
	path (refFasta)	

	output:
	path(${refFasta}.baseName.rpt_mask.gz), emit: refFasta.baseName

	script:
	"""
    genRefMask.py -r $refFasta -m 200 -p 95
    bgzip -c ${refFasta}.rpt.regions > ${refFasta.baseName}
	echo '##INFO=<ID=RPT,Number=1,Type=Integer,Description="Flag for variant in repetitive region">' > ${refFasta.baseName}.rpt_mask.hdr
	tabix -s1 -b2 -e3 ${refFasta.baseName}.rpt_mask.gz
    """
}

process BWA {
	cpus 8
	tag { "map clean ${uuid} reads to reference" }
    
	input:
	tuple val(uuid), path(reads), path(ref_index)

    output:
    tuple val(uuid), path("${uuid}.aligned.bam"), emit: bwa_mapped
    
	//don't add read group header here results in poorly formatted header
    
	"""
    bwa mem -r 1.5 -O 6 -t ${task.cpus} ref_index ${uuid}_clean.1.fq.gz ${uuid}_clean.2.fq.gz | samtools view -uS - | samtools sort -
	> bwa_mapped
    """
}

//return

//remove duplicates using samtools v 1.9
process removeDuplicates {
	tag "remove duplicates ${uuid}"
	
	publishDir "$params.outdir/${uuid}/bwa_mapped/${refFasta.baseName}/bam", mode: 'copy', pattern: "${uuid}.ba*"
    
	input:
    	tuple val (uuid), path (bwa_mapped)
    	path(ref_index)

    output:
    	tuple val(uuid), path("${uuid}.bam"), path("${uuid}.bam.bai") //emit: dup_removed

	//sort by name to run fixmate (to remove BWA artefacts) and provide info for markdup
	//sort by position to run markdup (and then remove duplicates)
    """
    samtools sort -@${task.cpus} -n -o sorted.bam ${uuid}.aligned.sam
    samtools fixmate -m sorted.bam fixed.sorted.bam
    samtools sort -@${task.cpus} -o fixed.resorted.bam fixed.sorted.bam
    samtools markdup -r fixed.resorted.bam ${uuid}.bam
    samtools index ${uuid}.bam
    """
}

//return

//run samtools mpileup - creates BCF containing genotype likelihoods 
process mpileup {

    input:
    	set uuid, file("${uuid}.bam"), file("${uuid}.bam.bai") 
    	file refFasta
    	file "*" from ref_index
 
    output:
    	set uuid, file("pileup.bcf") into pileup
   
    tag "${getShortId(uuid)}"
    //publishDir "$outputPath/$uuid/bwa_mapped/${refFasta.baseName}/vcf", mode: 'copy'

	//use bcftools mpileup to generate vcf file
	//mpileup genearates the likelihood of each base at each site
 	"""
    bcftools mpileup -Q25 -q30 -E -o40 -e20 -h100 -m2 -F0.002 -Ou -f $refFasta ${uuid}.bam > pileup.bcf
    """

}


//call SNPs using samtools call from mpileup file
process snpCall {

    input:
    	set uuid, file("pileup.bcf") from pileup
    	file refFasta
    	file "*" from ref_index
 
    output:
    	set uuid, file("${uuid}.bcf"), file("${uuid}.allsites.bcf") into snps_called
   
    tag "${getShortId(uuid)}"
    publishDir "$outdir/$uuid/bwa_mapped/${refFasta.baseName}/vcf", mode: 'copy'

	//call converts pileup to actual variants in the BCF or VCF file
	//norm, normalises indels
		//-m +any option to join biallelic sites into multiallelic records
		// and do this for any (i.e. SNPs and indels)
	
    """      
    # call variants only
    # 	-m use multiallelic model
    # 	-v output variants only
    bcftools call --prior 0.01 -Ou -m -v pileup.bcf | \
    	bcftools norm -f $refFasta -m +any -Ou -o ${uuid}.bcf
    	
    # call all sites
    bcftools call -Ou -m pileup.bcf | \
    	bcftools norm -f $refFasta -m +any -Ou -o ${uuid}.allsites.bcf
    """

}


//cleaner SNPS
process filterSnps {

    input:
    	set uuid, file("${uuid}.bcf"), file("${uuid}.allsites.bcf") 
    	file refFasta
    	file "*" from ref_index
 
    output:
    	set uuid, file("${uuid}.snps.vcf.gz"), file("${uuid}.snps.vcf.gz.csi"), 
    		file("${uuid}.zero_coverage.vcf.gz"), file("${uuid}.zero_coverage.vcf.gz.csi"), emit: filtered_snps
    	file "*"
   
    tag "${getShortId(uuid)}"
	publishDir "$outputPath/$uuid/bwa_mapped/${refFasta.baseName}/vcf", mode: 'copy', pattern: "${uuid}.snps.*"
	publishDir "$outputPath/$uuid/bwa_mapped/${refFasta.baseName}/vcf", mode: 'copy', pattern: "${uuid}.indels.*"
	publishDir "$outputPath/$uuid/bwa_mapped/${refFasta.baseName}/vcf", mode: 'copy', pattern: "${uuid}.zero_coverage.*"
	publishDir "$outputPath/$uuid/bwa_mapped/${refFasta.baseName}/vcf", mode: 'copy', pattern: "${uuid}.all.*"
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
	bcftools annotate -a ${refFasta.baseName}.rpt_mask.gz -c CHROM,FROM,TO,RPT 
		-h ${refFasta.baseName}.rpt_mask.hdr ${uuid}.bcf -Ob -o ${uuid}.masked.bcf.gz
    
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

//generate consensus fasta file
process consensusFa {

	input:
		set uuid, file("${uuid}.snps.vcf.gz"), file("${uuid}.snps.vcf.gz.csi"), 
    		file("${uuid}.zero_coverage.vcf.gz"), file("${uuid}.zero_coverage.vcf.gz.csi") from filtered_snps
		file refFasta
	
	output:
		set uuid, file("${uuid}.fa") into fa_file
		//file("*")
	
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
	cat tmp.fa | bcftools consensus -H 1 -M "-" ${uuid}.zero_coverage.vcf.gz > ${uuid}.fa
	"""
}


