#!/usr/bin/env python3

import argparse
import gzip
import pandas as pd
from Bio import SeqIO

def count_n_in_fasta(fasta_file, outfile):
    stats = {
        'seq_len' : 0,
        'N_count' : 0,
        'gap_count' : 0,
        'N_pc' : 0,
        'gap_pc' : 0,
        'pc_called' : 0,
    }

    with gzip.open(fasta_file, 'rt') if fasta_file.endswith('.gz') else open(fasta_file, 'r') as file:
        for record in SeqIO.parse(file, 'fasta'):
            sequence = str(record.seq)
            stats['seq_len'] += len(sequence)
            stats['N_count'] += sequence.count('N')
            stats['gap_count'] += sequence.count('-')

    stats['N_pc'] = round(100 * stats['N_count'] / stats['seq_len'], 2)
    stats['gap_pc'] = round(100 * stats['gap_count'] / stats['seq_len'], 2)
    stats['pc_called'] = round(100 - stats['N_pc'] - stats['gap_pc'], 2)

    df = pd.DataFrame([stats])
    df.to_csv(outfile, sep='\t', index=False)

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Count occurrences of "N" characters in a FASTA file.')
    parser.add_argument('fasta_file', help='Input FASTA file (can be gzipped)')
    parser.add_argument('outfile', help='Output file tsv')
    args = parser.parse_args()

    count_n_in_fasta(args.fasta_file, args.outfile)
