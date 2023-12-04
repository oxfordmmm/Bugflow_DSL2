
# This script will find the 7 gene loci in the reference genome
# Should be run from the directory containing the 7 gene loci
# Usage: bash find_loci_in_ref.sh <ref_genome> <output_file>

set -e

# get inputs
if [ -z "$1" ]
then
    echo "Please provide a reference genome"
    exit 1
fi

if [ -z "$2" ]
then
    echo "Please provide an output file"
    exit 1
fi


echo -e "chrom\tloci\tlength\tstart\tend" > $2
for gene in *.fas
do
    name=$(basename $gene .fas)
    echo $name

    blastn -query $gene -subject $1 \
        -out tmp.tsv -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore slen"

    echo -e "qseqid\tsseqid\tpident\tlength\tmismatch\tgapopen\tqstart\tqend\tsstart\tsend\tevalue\tbitscore\tslen" > $name.tsv && cat tmp.tsv >> $name.tsv

    rm tmp.tsv

    python3 get_best_match.py -t $name.tsv >> $2

    rm $name.tsv
done




