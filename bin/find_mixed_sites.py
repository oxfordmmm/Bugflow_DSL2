#!/usr/bin/env python3

import argparse
import pandas as pd
import pysam

# Indexing:
# blastn outputs are 1-indexed and ranges are inclusive
# i.e. adk_1 is 501 bases and goes from 118604 to 119104 inclusive
# Pysam pileup is 0-indexed though!



# Set the mpileup options
mpileup_options = {
    "min_base_quality": 30,
    "min_mapping_quality": 30,
    "redo_baq": True
}

def Bases_At_Pos(samfile, pos, chromname):
    'Return a string of the bases at that position (1-indexed).'
    position = 0
    coverage = 0
    bases = ""

    # Need to convert request to o-index
    for pileupcolumn in samfile.pileup(reference=chromname, start=pos-1, end=pos, truncate=True, **mpileup_options):
        position = int(pileupcolumn.pos+1)  # converting back to 1-based
        coverage = int(pileupcolumn.n)
        for pileupread in pileupcolumn.pileups:
            if (pileupread.indel == 0 and pileupread.is_del == 0):
                bases += pileupread.alignment.seq[pileupread.query_position]
    return position, coverage, bases

def eyre_method(bam, outfile, loci_df: pd.DataFrame):
    samfile = pysam.Samfile(bam, 'rb')

    with open (outfile, 'w') as f:
        f.write('site\tA\tC\tG\tT\n')
        #for each position of interest write base counts
        for chrom, start, end in zip(loci_df['chrom'], loci_df['start'], loci_df['end']):
            for pos in range(start, end+1):
                hq_position, _, hq_bases = Bases_At_Pos(samfile, int(pos), chrom)
                hq_basecounts = [ len([x for x in hq_bases if x == B]) for B in ['A', 'C', 'G', 'T'] ]
                if sum(hq_basecounts) > 0:
                    f.write('%s\t%s\n'%(hq_position, '\t'.join([str(i) for i in hq_basecounts])))

        f.close()


def pysam_pileup(bam, ref, chrom, loci, start, end):
    # convert request to 0-base
    start = start-1
    end = end-1

    # Perform mpileup from start to end inclusive
    pileup = bam.pileup(chrom, start, end+1, truncate=True, **mpileup_options)
    variant_sites = []

    for pileupcolumn in pileup:
        reads_all = len(pileupcolumn.pileups)
        # note: pileupcolumn.n ignores the quality filter so returns a higher number than expected

        column = {
            'loci' : loci,
            'chrom' : chrom,
            'pos' : pileupcolumn.pos + 1, # return to 1-based indexing
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
        
        
        column['loci'] = loci
        column['major'] = major
        column['major_pc'] = round(100 * major / reads_all, 1) if reads_all > 0 else 0.0
        column['rest'] = rest
        column['rest_pc'] = round(100 * rest / reads_all, 1) if reads_all > 0 else 0.0
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
    parser.add_argument("--loci", required=True, help="Path to the TSV file containing chrom, start, and end positions of relevant loci")
    parser.add_argument("--outfile_prefix", required=True, help="Prefix for output csv files")
    args = parser.parse_args()

    ref = pysam.FastaFile(args.ref)
    bam = pysam.AlignmentFile(args.bam)
    loci_df = pd.read_csv(args.loci, sep='\t', dtype={'chrom': str, 'loci': str, 'start': int, 'end': int})
    outfile_prefix = args.outfile_prefix

    eyre_method(args.bam, f"{outfile_prefix}_eyre.tsv", loci_df)

    variant_sites = []
    total_sites = 0
    for chrom, loci, start, end in zip(loci_df['chrom'], loci_df['loci'], loci_df['start'], loci_df['end']):
        # For each loci the function returns a list of dictionaries with a dict for each column in the pileup
        sites = pysam_pileup(bam, ref, chrom, loci, start, end)
        variant_sites.extend(sites)
        total_sites += end-start
    
    col_names = ['loci', 'chrom', 'pos', 'reads_all', 'major', 'rest', 'major_pc', 'rest_pc', 'insertions', 'deletions', 'A', 'C', 'T', 'G']
    df = pd.DataFrame(columns=col_names)
    if len(variant_sites) > 0:
        df = pd.DataFrame(variant_sites)[col_names]
    df.to_csv(f"{outfile_prefix}.csv", index=False)
    
    summary_df = summarise_mixed_sites(df, total_sites)
    summary_df.to_csv(f"{outfile_prefix}_summarised.csv", index=False)

if __name__ == "__main__":
    main()
