#!/usr/bin/env python3

import argparse
import pandas as pd
from Bio import SeqIO
from Bio.Seq import Seq, MutableSeq
import itertools

def add_extra_characters(fasta_file, snps_vcf, allsites_vcf, repeat_mask, output_file):
    # Read the fasta file
    fasta_data = SeqIO.read(fasta_file, 'fasta')
    seq = MutableSeq(fasta_data.seq)
    ref_length = len(seq)
    print(f"Reference length: {ref_length}")

    # Find all zero depth positions, according to allsites vcf
    allsites_df = pd.read_csv(allsites_vcf, sep='\t', comment="#", header=None,
                           names=['chrome', 'POS', 'id', 'ref', 'alt', 'qual', 'filter', 'info', 'format', 'sample'],
                           dtype=str)
    allsites_df['POS'] = allsites_df['POS'].astype(int)

    # Filter out rows where 'info' contains "DP=0"
    allsites_df = allsites_df[~allsites_df['info'].str.contains("DP=0")]

    nonzero_positions = set(allsites_df['POS'])
    zero_positions = set(range(1, ref_length + 1)) - nonzero_positions

    # Get repeat regions from mask
    repeat_df = pd.read_csv(repeat_mask, sep='\t', comment="#", header=None,
                            names=['chrome', 'start', 'end', 'RPT'])
    ranges = [(start, end) for start, end in zip(repeat_df['start'], repeat_df['end'])]
    repeat_positions = set(itertools.chain.from_iterable(
        [range(start, end+1) for start, end in ranges]
    ))

    snps_df = pd.read_csv(snps_vcf, sep='\t', comment="#", header=None,
                          names=['chrome', 'POS', 'id', 'ref', 'alt', 'qual', 
                                 'filter', 'info', 'format', 'sample'])
    filtered_locations = set(snps_df[snps_df['filter'] != 'PASS']['POS'])
    het_positions = set(snps_df[snps_df['filter'].str.contains('HetroZ')]['POS'])
    null_gt_positions = set(snps_df[snps_df['sample'].str.contains('./.', regex=False)]['POS'])


    # Replace in order of precedence
    already_changed = set()
    position_lists = [repeat_positions, zero_positions, null_gt_positions, het_positions, filtered_locations]
    characters = ['R', 'Z', 'X', 'H', 'F']
    descriptions = ['repeat', 'zero depth', 'null genotype', 'heterozygous', 'filtered']
    for pos_list, character, desc in zip(position_lists, characters, descriptions):
        pos_list -= already_changed
        print(f"{desc} ({character}): {len(pos_list)} ({100 * len(pos_list)/ref_length:.2f}%)")
        for position in pos_list:
            seq[position-1] = character
        already_changed = already_changed | pos_list
    
    fasta_data.seq = Seq(str(seq))

    # Write the modified fasta sequence to a new file
    with open(output_file, 'w') as file:
        SeqIO.write(fasta_data, file, 'fasta')

if __name__ == '__main__':
    # Parse command-line arguments
    parser = argparse.ArgumentParser(description='Mask positions with depth 0 in a fasta file')
    parser.add_argument('-f', '--fasta_file', help='Path to the fasta file')
    parser.add_argument('-a', '--allsites_vcf', help='Path to the allsites vcf file')
    parser.add_argument('-s', '--snps_vcf', help='Path to the snps vcf file')
    parser.add_argument('-r', '--repeat_mask', help='Path to the repeat mask file')
    parser.add_argument('-o', '--output_file', help='Path to the output file')
    args = parser.parse_args()

    # Call the mask_zero_depth function
    add_extra_characters(args.fasta_file, args.snps_vcf, args.allsites_vcf,
                    args.repeat_mask, args.output_file)
