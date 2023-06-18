#!/usr/bin/env python3

import argparse
import pandas as pd
import pysam
import pysamstats

def pysam_pileup(bam, ref, chrom, loci, start, end):

    # Set the mpileup options
    mpileup_options = {
        "min_base_quality": 25,
        "min_mapping_quality": 30,
        "redo_baq": True
    }

    # Perform mpileup
    # pileup = bam.pileup(bam, **mpileup_options)
    pileup = bam.pileup(chrom, start, end, truncate=True, **mpileup_options)
    variant_sites = []

    for pileupcolumn in pileup:
        # print("\ncoverage at base %s = %s" % (pileupcolumn.pos, pileupcolumn.n))
        reads_all = len(pileupcolumn.pileups)
        # note: pileupcolumn.n ignores the quality filter so returns a higher number than expected

        column = {
            'loci' : loci,
            'chrom' : chrom,
            'pos' : pileupcolumn.pos,
            'reads_all' : reads_all,
            'insertions' : 0,
            'deletions' : 0,
            'A' : 0,
            'C' : 0,
            'T' : 0,
            'G' : 0,
        }

        for pileupread in pileupcolumn.pileups:
            if pileupread.indel > 0:
                column['insertions'] += 1
                continue
            if pileupread.is_del > 0:
                column['deletions'] += 1
                continue

            base = pileupread.alignment.query_sequence[pileupread.query_position]
            column[base] += 1
        
        # check if mixed
        calls = [column['insertions'], column['deletions'], column['A'], column['C'], column['T'], column['G']]
        sorted_count = sorted(calls, reverse=True)
        
        major = sorted_count[0]
        rest = sum(sorted_count[1:])
        
        if rest > 2 or column['insertions'] > 0 or column['deletions'] > 0:
            column['loci'] = loci
            column['major'] = major
            column['major_pc'] = round(100 * major / reads_all, 1)
            column['rest'] = rest
            column['rest_pc'] = round(100 * rest / reads_all, 1)
            variant_sites.append(column)

    return variant_sites


# Simpler version using pysamstats, but can't filter the reads based on quality
def find_variant_sites(bam, fasta, chromosome, loci, start, end):
    stats_variation = pysamstats.stat_variation(bam, fasta, chrom=chromosome, start=start, end=end, truncate=True)
    
    variant_sites = []
    for _column in stats_variation:
        column = _column.copy()
    
        reads_all = column['reads_all']
        insertions = column['insertions']
        deletions = column['deletions']

        # Check if there are more than two bases in the minor group
        calls = [insertions, deletions, column['A'], column['C'], column['T'], column['G']]
        sorted_count = sorted(calls, reverse=True)
        
        major = sorted_count[0]
        rest = sum(sorted_count[1:])
        
        if rest > 2 or insertions > 0 or deletions > 0:
            column['loci'] = loci
            column['major'] = major
            column['major_pc'] = round(100 * major / reads_all, 1)
            column['rest'] = rest
            column['rest_pc'] = round(100 * rest / reads_all, 1)
            variant_sites.append(column)
    return variant_sites

def summarise_mixed_sites(df, total_bases):
    # want total number of sites, number 
    rest_1_count = (df['rest'] == 1).sum()
    rest_n_count = (df['rest'] >= 2).sum()
    other_base_1_count = (df['rest'] - df['insertions'] - df['deletions'] == 1).sum()
    other_base_n_count = (df['rest'] - df['insertions'] - df['deletions'] >= 2).sum()
    del_1_count = (df['deletions'] == 1).sum()
    del_n_count = (df['deletions'] >= 2).sum()
    insert_1_count = (df['insertions'] == 1).sum()
    insert_n_count = (df['insertions'] >= 2).sum()

    values = {
        'total_positions': total_bases,
        'rest_1' : rest_1_count,
        'rest_2_plus' : rest_n_count,
        'other_base_1' : other_base_1_count,
        'other_base_2_plus' : other_base_n_count,
        'del_1' : del_1_count,
        'del_2_plus' : del_n_count,
        'insert_1' : insert_1_count,
        'insert_2_plus' : insert_n_count,
    }
    
    df = pd.DataFrame([values])
    return df



def main():
    parser = argparse.ArgumentParser(description="Calculate variation statistics and check for positions with a mix of bases called.")
    parser.add_argument("--ref", required=True, help="Path to the reference FASTA file")
    parser.add_argument("--bam", required=True, help="Path to the BAM file")
    parser.add_argument("--loci", required=True, help="Path to the TSV file containing start and end positions of relevant loci")
    parser.add_argument("--outfile_prefix", required=True, help="Prefix for output csv files")
    parser.add_argument("--chrom", required=False, default='CP010905.2', help="Chromosome on ref FASTA file")
    args = parser.parse_args()

    ref = pysam.FastaFile(args.ref)
    bam = pysam.AlignmentFile(args.bam)
    loci_file = args.loci
    chrom = args.chrom
    outfile_prefix = args.outfile_prefix

    loci_df = pd.read_csv(loci_file, sep='\t')
    variant_sites = []
    total_sites = 0
    for loci, start, end in zip(loci_df['loci'], loci_df['start'], loci_df['end']):
        # For each loci the function returns a list of dictionaries with a dict for each column in the pileup
        sites = pysam_pileup(bam, ref, chrom, loci, int(start), int(end))
        variant_sites.extend(sites)
        total_sites += end-start
    
    col_names = ['loci', 'chrom', 'pos', 'reads_all', 'major', 'rest', 'major_pc', 'rest_pc', 'insertions', 'deletions', 'A', 'C', 'T', 'G']
    df = pd.DataFrame(columns=col_names)
    if len(variant_sites) > 0:
        df = pd.DataFrame(variant_sites)
        df = df[col_names]
    df.to_csv(f"{outfile_prefix}.csv", index=False)
    
    summary_df = summarise_mixed_sites(df, total_sites)
    summary_df.to_csv(f"{outfile_prefix}_summarised.csv", index=False)

if __name__ == "__main__":
    main()
