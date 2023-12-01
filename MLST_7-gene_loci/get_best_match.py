import argparse
import pandas as pd

if __name__ == '__main__':
    # Parse command-line arguments
    parser = argparse.ArgumentParser(description='Get the best match from a TSV file')
    parser.add_argument('-t', '--tsv', help='Path to the TSV file')
    args = parser.parse_args()
    
    # Read the TSV file using pandas
    df = pd.read_csv(args.tsv, sep='\t')

    # Find the row with the highest 'bitscore'
    max_bitscore_row = df.loc[df['bitscore'].idxmax()]

    start = min(max_bitscore_row['sstart'], max_bitscore_row['send'])
    end = max(max_bitscore_row['sstart'], max_bitscore_row['send'])

    print(f"{max_bitscore_row['qseqid']}\t{max_bitscore_row['length']}\t{start}\t{end}")