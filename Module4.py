# SVModeller - Module 4

# Modify reference genome

# Input:
# - Events to mofidy reference genome (.tsv)
# - Reference genome (.fasta)
# - OPTIONAL: additional events table (.tsv)

# Output:
# - Modified reference genome (Modified_Reference_Genome.fasta)
# - Table with added events and their current positions in the genome (Sorted_Genomic_Events.tsv)

# Developers
# SVModeller has been developed by Ismael Vera-Munoz (orcid.org/0009-0009-2860-378X) at the Repetitive DNA Biology (REPBIO) Lab at the Centre for Genomic Regulation (CRG) (Barcelona 2024-2026)

# License
# SVModeller is distributed under the AGPL-3.0.

import argparse
import pandas as pd
from GAPI import formats
import warnings

def write_fasta(file_path, seq_dict):
    '''Create a FASTA file from a dictionary of sequences'''
    with open(file_path, 'w') as fasta_file:
        for header, sequence in seq_dict.items():
            fasta_file.write(f">{header}\n")
            fasta_file.write(f"{sequence}\n")

# Remove FutureWarnings
warnings.simplefilter(action='ignore', category=FutureWarning)

def main(file1, file2, fasta_file):
    print(f'File 1 path: {file1}')
    if file2:
        print(f'File 2 path: {file2}')
    print(f'FASTA file: {fasta_file}')

    # Load and merge input TSV files
    df_1 = pd.read_csv(file1, sep='\t')
    merged_df = pd.concat([df_1, pd.read_csv(file2, sep='\t')], ignore_index=True) if file2 else df_1
    sorted_df = merged_df.sort_values(by=['#ref', 'beg', 'Event_Type'])

    # Read reference genome
    fasta_reader = formats.FASTA()
    fasta_reader.read(fasta_file)

    # Initialize offsets and columns for haplotype positions
    offsets = {ref: 0 for ref in fasta_reader.seqDict.keys()}
    sorted_df.insert(4, 'beg_haplotype', None)
    sorted_df.insert(5, 'end_haplotype', None)

    # Apply events and calculate haplotype coordinates on the fly
    for index, row in sorted_df.iterrows():
        ref_key = row['#ref'].strip()
        beg = int(row['beg'])
        event = row['Event_Type']
        length = int(row['Length']) if pd.notna(row['Length']) else 0
        seq_insertion = row['Sequence_Insertion'] if 'Sequence_Insertion' in row and pd.notna(row['Sequence_Insertion']) else ''

        if ref_key not in fasta_reader.seqDict:
            continue

        # Adjust position based on current offset
        adjusted_beg = beg + offsets[ref_key]

        if event == 'Insertion':
            # Modify reference coordinates: end = beg + 1
            sorted_df.at[index, 'end'] = beg + 1  # Insertions should have end = beg + 1
            
            # Do not modify the Length here; keep the original length of the insertion
            # sorted_df.at[index, 'Length'] = 1  # Removed this line to keep the original length

            # Insert into sequence
            if adjusted_beg <= len(fasta_reader.seqDict[ref_key]):
                fasta_reader.seqDict[ref_key] = (
                    fasta_reader.seqDict[ref_key][:adjusted_beg] +
                    seq_insertion +
                    fasta_reader.seqDict[ref_key][adjusted_beg:]
                )
                # Update offset based on the length of the inserted sequence
                offsets[ref_key] += len(seq_insertion)

            # Haplotype coordinates for insertion
            sorted_df.at[index, 'beg_haplotype'] = adjusted_beg
            sorted_df.at[index, 'end_haplotype'] = adjusted_beg + len(seq_insertion)

        elif event == 'Deletion':
            if adjusted_beg < len(fasta_reader.seqDict[ref_key]):
                fasta_reader.seqDict[ref_key] = (
                    fasta_reader.seqDict[ref_key][:adjusted_beg] +
                    fasta_reader.seqDict[ref_key][adjusted_beg + length:]
                )
                # Update offset based on the deletion length
                offsets[ref_key] -= length

            # Haplotype coordinates for deletion (only 1bp apart)
            sorted_df.at[index, 'beg_haplotype'] = adjusted_beg
            sorted_df.at[index, 'end_haplotype'] = adjusted_beg + 1  # Only 1bp apart for deletions

    # Save outputs
    sorted_df.to_csv('Sorted_Genomic_Events.tsv', sep='\t', index=False)
    write_fasta("Modified_Reference_Genome.fasta", fasta_reader.seqDict)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Process genomic insertions and deletions')
    parser.add_argument('--file1', type=str, required=True, help='Path to the first TSV file.')
    parser.add_argument('--fasta_file', type=str, required=True, help='Path to the FASTA file.')
    parser.add_argument('--file2', type=str, nargs='?', default=None, required=False, help='Optional path to the second TSV file.')
    args = parser.parse_args()
    main(args.file1, args.file2, args.fasta_file)
