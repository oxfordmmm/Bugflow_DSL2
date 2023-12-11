# Parse BLASTN output and print the best match for each gene
import pandas as pd
import argparse
from Bio import SeqIO

def main():
    # Get arguements using argparse
    parser = argparse.ArgumentParser(description=
        'Parse BLASTN output and print the best match for each gene')
    parser.add_argument('--assembly', help='Assembly fasta')
    parser.add_argument('--blastn_output', help='BLASTN output file')
    parser.add_argument('--output', help='Output file name')
    args = parser.parse_args()

    # Read BLASTN output into a pandas DataFrame
    blast_results = pd.read_csv(args.blastn_output, sep='\t')

    # Group by query and find the best match for each gene based on the highest bit score
    best_matches = blast_results.groupby('qseqid').apply(lambda x: x.loc[x['bitscore'].idxmax()])
    best_matches['start'] = best_matches[['sstart', 'send']].min(axis=1)
    best_matches['end'] = best_matches[['sstart', 'send']].max(axis=1)
    best_matches['reverse'] = best_matches['sstart'] > best_matches['send']
    
    # Print the results
    print(best_matches)

    assembly = SeqIO.parse(open(args.assembly, 'r'), 'fasta')
    genome_dict = SeqIO.to_dict(assembly)

    with open(args.output, 'w') as output_file:
        for gene, contig, start, end, reverse in zip(best_matches['qseqid'], best_matches['sseqid'], best_matches['start'], best_matches['end'], best_matches['reverse']):
            # blast is 1-indexed so adjust
            start -= 1
            end -= 1
            print(gene, contig, start, end, reverse)
            match = genome_dict[contig].seq[start:end + 1]
            if reverse:
                match = match.reverse_complement()
            
            print(match[:10])
            output_file.write(f'>{gene}\n{match}\n')

if __name__ == '__main__':
    main()