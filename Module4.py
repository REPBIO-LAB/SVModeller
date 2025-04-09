# SVModeller - Module 4

# Modify reference genome

# Input:
# - Events to mofidy reference genome (.tsv)
# - Reference genome (.fasta)
# - OPTIONAL: additional events table (.tsv)

# Output:
# - Modified reference genome (Modified_Reference_Genome.fasta)
# - Table with added events and their current positions in the genome (Sorted_Genomic_Events.tsv)

import argparse
import pandas as pd
import formats
import os

def write_fasta(file_path, seq_dict):
    ''' 
    Function to create a fasta file from a dictionary
    '''
    with open(file_path, 'w') as fasta_file:
        for header, sequence in seq_dict.items():
            # Write the header line
            fasta_file.write(f">{header}\n")  
            # Write the sequence line
            fasta_file.write(f"{sequence}\n")  

def main(file1, file2, fasta_file):
    # Print the paths of the input files
    print(f'File 1 path: {file1}')
    if file2:
        print(f'File 2 path: {file2}')
    print(f'FASTA file: {fasta_file}')

    # Load the first file (insertions file)
    df_1 = pd.read_csv(file1, sep='\t')
    
    # If a second file is provided, load and merge
    if file2:
        df_2 = pd.read_csv(file2, sep='\t')
        # Merging the DataFrames by concatenating them vertically
        merged_df = pd.concat([df_1, df_2], ignore_index=True, sort=False)
    else:
        # If no second file, just use insertions
        merged_df = df_1

    # Order according to chromosome, start point, and event type
    sorted_df = merged_df.sort_values(by=['#ref', 'beg', 'Event_Type'])

    # Create an instance of the FASTA class
    fasta_reader = formats.FASTA()
    # Read the FASTA file
    fasta_reader.read(fasta_file)

    # Dictionary with keys 'chr1' to 'chr23' and all values set to 0
    counter_dict = {f'chr{i}': 0 for i in range(1, 23)}
    counter_dict['chrX'] = 0
    counter_dict['chrY'] = 0

    # For each row get the necessary values
    for index, row in sorted_df.iterrows():
        ref_key = row['#ref']
        beg_position = row['beg']
        total_length = row['Length']
        event_type = row['Event_Type']
        
        if ref_key in fasta_reader.seqDict:
            # Insertions
            if event_type == 'Insertion':
                # Get the sequence
                sequence_ins = row['Sequence_Insertion']
                if beg_position <= len(fasta_reader.seqDict[ref_key]):
                    # Divide the sequence in 2 and insert the sequence
                    fasta_reader.seqDict[ref_key] = (
                        fasta_reader.seqDict[ref_key][:beg_position] +
                        sequence_ins +
                        fasta_reader.seqDict[ref_key][beg_position:]
                    )
                # Sum the length to the counter
                counter_dict[ref_key] += total_length
            
            # Deletions
            elif event_type == 'Deletion':
                if beg_position < len(fasta_reader.seqDict[ref_key]):
                    # Split the sequence in 2, removing the part of the deletion
                    fasta_reader.seqDict[ref_key] = (
                        fasta_reader.seqDict[ref_key][:beg_position] +
                        fasta_reader.seqDict[ref_key][beg_position + total_length:]
                    )
                # Subtract the value from the counter
                counter_dict[ref_key] -= total_length

    # Initialize the new columns in the DataFrame
    sorted_df.insert(4, 'beg_haplotype', None)
    sorted_df.insert(5, 'end_haplotype', None)

    # Maintain a dynamic starting length based on the counter before changes
    for index, row in sorted_df.iterrows():
        ref_key = row['#ref']
        beg_position = row['beg']
        total_length = row['Length']
        event_type = row['Event_Type']
        
        # Calculate current sequence length based on counter
        current_length = counter_dict[ref_key]
        
        # Insertions
        if event_type == 'Insertion':
            # Get the sequence
            sequence_ins = row['Sequence_Insertion']
            # Divide the sequence in 2 and insert the sequence
            fasta_reader.seqDict[ref_key] = (
                fasta_reader.seqDict[ref_key][:beg_position] +
                sequence_ins +
                fasta_reader.seqDict[ref_key][beg_position:]
            )
            # Update the new columns
            sorted_df.at[index, 'beg_haplotype'] = (
                current_length + len(fasta_reader.seqDict[ref_key][:beg_position])
            )
            sorted_df.at[index, 'end_haplotype'] = (
                sorted_df.at[index, 'beg_haplotype'] + len(sequence_ins)
            )

            # Increase the length of the sequence in the counter
            counter_dict[ref_key] += total_length
            
        # Deletions
        elif event_type == 'Deletion':
            # Calculate start and end positions for deletion
            start_position = beg_position
            end_position = beg_position + total_length 
            
            # Update the new columns 
            start = current_length + len(fasta_reader.seqDict[ref_key][:beg_position])
            sorted_df.at[index, 'beg_haplotype'] = start
            sorted_df.at[index, 'end_haplotype'] = start + 1
            
            # Split the sequence in 2, removing the part of the deletion
            fasta_reader.seqDict[ref_key] = (
                fasta_reader.seqDict[ref_key][:start_position] +
                fasta_reader.seqDict[ref_key][end_position:]
            )

            # Decrease the length of the sequence in the counter
            counter_dict[ref_key] -= total_length
    
    # Save table with positions after modifications
    sorted_df.to_csv('Sorted_Genomic_Events.tsv', sep='\t', index=False)
    # Write the modified sequences to a new FASTA file
    write_fasta("Modified_Reference_Genome.fasta", fasta_reader.seqDict)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Process genomic insertions and deletions')
    parser.add_argument('file1', type=str, help='Path to the first TSV file .')
    parser.add_argument('fasta_file', type=str, help='Path to the FASTA file.')
    parser.add_argument('file2', type=str, nargs='?', default=None, help='Optional path to the second TSV file.')
    args = parser.parse_args()
    main(args.file1, args.file2, args.fasta_file)
