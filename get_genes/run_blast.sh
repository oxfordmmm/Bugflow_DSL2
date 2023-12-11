
for assembly in consensus_fa/*.fa
do
    name=$(basename $assembly .fa)
    echo $name

    blastn -query Cdiff_AMR.fasta -subject $assembly \
        -out $name.tmp.tsv -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore slen" \
        -subject_besthit

    echo -e "qseqid\tsseqid\tpident\tlength\tmismatch\tgapopen\tqstart\tqend\tsstart\tsend\tevalue\tbitscore\tslen" > $name.blastn.tsv && cat $name.tmp.tsv >> $name.blastn.tsv

    python3 get_best.py --assembly $assembly --blastn_output $name.blastn.tsv --output $name.matches.fa
done



