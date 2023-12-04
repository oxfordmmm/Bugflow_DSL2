source /mnt/scratch/miniconda3/etc/profile.d/conda.sh

batch=batch_name
echo running $batch

input='<path to input>/*_R{1,2}.fastq.gz'
ref="refs/Cdiff_630_GCA_000009205.1.fasta"
ref_mlst="refs/Cdiff_630_GCA_000009205.1.mlst_loci.tsv"
run_dir=runs/$batch
output_dir=output/$batch


mkdir -p $run_dir/cgmlst $run_dir/mapping

cd $run_dir/cgmlst
nextflow run Bugflow_DSL2 -entry cdiff_hcgmlst_amrg_blastn_single --reads $input \
    --outdir $output_dir/cgmlst -profile docker,oci

cd $run_dir/mapping
nextflow run Bugflow_DSL2 -entry cdiff_mapping_snpCalling_DE --reads $input \
    --refFasta $ref --mlst_loci $ref_mlst --outdir $output_dir/mapping -profile docker,oci
