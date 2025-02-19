# SVMoldeler - Module 3

# Process Deletions

import argparse
import formats
import gRanges
import structures
import numpy as np
import pandas as pd

def read_vcf_file_BED(file_path):
    ''' 
    Funtion to:
    - Inicializate a VCF object
    - use the VCF class to read the path to the VCF file
    - Create empty list to store the data
    - For each line get the information and place it in a new dataframe
    - First the values that are not common in all lines are defined to get their value or write NA
    - Add to the list all the information and add this list to the dataframe 
    '''
    data = [] # List to store the extracted data
    with open(file_path, 'r') as f:
        for line in f:
            if line.startswith('#'):
                continue
            else:
                elements = line.strip().split('\t')
                
                # Get the identifier from the third element
                identifier = elements[2]
                
                # Get the chromosome and position from the first and second elements
                chromosome = elements[0]
                position = int(elements[1])
                
                # Get the info field
                info = elements[7].split(';')

                family, strand, length, dtype_n = ('NA',) * 4
                for item in info:
                    key, value = item.split('=')
                    if key == 'FAM_N':
                        family = value
                    elif key == 'STRAND':
                        strand = value
                    elif key == 'INS_LEN':
                        length = value
                    elif key == 'DTYPE_N':
                        dtype_n = value

                end_position = position + int(length) if length != 'NA' else 'NA'
                seq_len = len(elements[3])

                # Add to the data list the information
                data.append({
                    'Identifier': identifier,
                    'Type_SV': dtype_n, 
                    'Family': family,
                    'Start_position': position,
                    'End_position': end_position,
                    'Length': seq_len,
                    'Chromosome': chromosome,
                    'Strand': strand
                })

    df = pd.DataFrame(data, columns=['Chromosome', 'Start_position', 'End_position', 'Identifier', 'Type_SV', 'Family', 'Length', 'Strand'])
    return df

def merge_dataframes(df1, df2):
    ''' 
    Function to merge two dataframes
    '''
    # Make a copy of the first DataFrame
    df1_copy = df1.copy()

    # Identify the rows in df2 where the 'identifier' is not in df1
    df2_filtered = df2[~df2['Identifier'].isin(df1['Identifier'])]

    # Append the filtered rows from df2 to the copy of df1
    return pd.concat([df1_copy, df2_filtered], ignore_index=True)

def mut_bins(bins, result_BED_table):
    """
    Classify mutations in each window of the bindatabase.

    Parameters:
    bins (list): List of windows in the format [(chromosome, start, end), ...]
    result_BED_table (pandas.DataFrame): Table with mutation data

    Returns:
    pandas.DataFrame: Classified mutations in each window
    """
    # Create a dictionary to map SubType to column names
    subtype_cols = {
        'solo_SVA': 'solo_SVA',
        'solo_Alu': 'solo_Alu',
        'solo_L1': 'solo_L1',
        'partnered_SVA': 'partnered_SVA',
        'partnered_L1': 'partnered_L1',
        'orphan': 'orphan',
        'simple': 'simple',
        'VNTR': 'VNTR'
    }

    # Initialize the output df
    df = pd.DataFrame(columns=['window', 'beg', 'end'] + list(subtype_cols.values()))

    # Iterate over each window
    for window in bins:
        chrom, start, end = window
        window_df = result_BED_table[(result_BED_table['#ref'] == chrom) & (result_BED_table['beg'] >= start) & (result_BED_table['end'] <= end)]

        # Count the occurrences of each SubType in the window
        counts = window_df['name'].value_counts().to_dict()

        # Create a row for the output DataFrame
        row = [chrom, start, end] + [counts.get(subtype, 0) for subtype in subtype_cols]

        # Append the row to the output df
        df.loc[len(df)] = row
    return df

def calculate_normalized_values(res_table):
    ''' 
    Function to calculate the normalized values of each variant of each chromosome
    '''
    # Get unique values in the 'window' column
    unique_chr_values = res_table['window'].unique()

    #new df to store the results
    new_table = pd.DataFrame(columns=res_table.columns)
    for chr_value in unique_chr_values:

        # Filter the rows for the current 'window' value
        subset = res_table[res_table['window'] == chr_value]

        # Copy the values of 'beg' and 'end' columns
        normalized_values = subset[['beg', 'end']]

        # Calculate the sum of values for each column
        total_values = subset[['solo_SVA', 'solo_Alu', 'solo_L1', 'partnered_SVA',
                                'partnered_L1', 'orphan', 'simple', 'VNTR']].sum()
        
        # Divide each cell value by the total sum
        normalized_values[['solo_SVA', 'solo_Alu', 'solo_L1', 'partnered_SVA',
                           'partnered_L1', 'orphan', 'simple', 'VNTR']] = subset[['solo_SVA', 'solo_Alu', 'solo_L1', 'partnered_SVA',
                                                                                    'partnered_L1', 'orphan', 'simple', 'VNTR']].div(total_values).values
        
        # Add the 'window' values back to the df
        normalized_values['window'] = chr_value

        # Append the normalized values to the new df
        new_table = pd.concat([new_table, normalized_values], ignore_index=True)
    return new_table

def normalize_columns(df):
    ''' 
    Function to normalize columns of a given dataframe
    '''
    columns_to_normalize = df.columns[3:]
    for column in columns_to_normalize:
        total_sum = df[column].sum()
        if total_sum != 0:
            df[column] = df[column] / total_sum
    return df

def add_columns(df1, df2, df3):
    ''' 
    Function to add start and end columns to the df based on probabilities
    '''
    new_df = pd.DataFrame(columns=['#ref', 'beg', 'end', 'Length', 'Strand'])
    for _, row in df3.iterrows():
        name = row['name']
        prob_df = df2[df2[name] > 0]
        if prob_df.empty:
            random_row = df1[df1['name'] == name].sample(n=1)
        else:
            name_df = df1[df1['name'] == name]
            if name_df.empty:
                continue
            weights = prob_df[name].values
            if weights.sum() == 0:
                weights = np.ones(len(weights)) / len(weights)
            weights = np.repeat(weights, len(name_df) // len(weights) + 1)[:len(name_df)]
            random_row = name_df.sample(n=1, weights=weights)
        new_df = new_df._append(random_row[['#ref', 'beg', 'end', 'Length', 'Strand']], ignore_index=True)
    df3[['#ref', 'beg', 'end', 'Length', 'Strand']] = new_df
    return df3

def main(vcf_path_all, vcf_path_class, bin_size, chromosome_length, num_events):
    # Print the paths of the input files
    print(f'VCF file with all deletions: {vcf_path_all}')
    print(f'VCF file with classified deletions: {vcf_path_class}')
    print(f'Chromosome length file: {chromosome_length}')
    print(f'Defined number of events: {num_events}')
    
    # Process VCF file with ALL deletions
    result_df = read_vcf_file_BED(vcf_path_all)
    result_BED_table_all_deletion = result_df.rename(columns={'Type_SV': 'name', 'Chromosome':'#ref', 'Start_position': 'beg', 'End_position':'end', 'Family':'SubType'})
    result_BED_table_all_deletion.fillna('NA', inplace=True)

    # Process VCF with classified deletions
    result_df = read_vcf_file_BED(vcf_path_class)
    result_BED_table_class_deletion = result_df.rename(columns={'Type_SV': 'name', 'Chromosome':'#ref', 'Start_position': 'beg', 'End_position':'end', 'Family':'SubType'})
    result_BED_table_class_deletion.fillna('NA', inplace=True)

    # Merge the two VCF files (processed)
    merged_df = merge_dataframes(result_BED_table_class_deletion, result_BED_table_all_deletion)

    # Add end column
    na_rows = merged_df['end'] == 'NA'
    merged_df.loc[na_rows, 'end'] = merged_df.loc[na_rows, 'beg'] + merged_df.loc[na_rows, 'Length']
    merged_df['end'] = merged_df['end'].astype(int)

    # Update the names of the events
    merged_df['name'] = merged_df['name'] + '_' + merged_df['SubType'].apply(lambda x: x if x in ['SVA', 'Alu', 'L1'] else '')
    merged_df = merged_df.drop('SubType', axis=1)
    merged_df['name'] = merged_df['name'].apply(lambda x: x[:-1] if x.endswith('_') else x)
    merged_df['name'] = merged_df['name'].replace('NA', 'simple')

    # Remove non desired ones
    merged_df = merged_df[~merged_df['name'].isin(['partnered_Alu', 'chimera', 'PSD', 'partnered_SVA', 'partnered_L1'])]

    # Chromosomes length
    chr_length = formats.chrom_lengths_index(chromosome_length)

    # Bins to fragment the genome
    bins = gRanges.makeGenomicBins(chr_length, bin_size, [f'chr{i}' for i in range(1, 23)])[::-1]

    # Classify mutations in each window of the bindatabase
    res_table = mut_bins(bins, merged_df)
    
    # Normalize values
    prob_INS = calculate_normalized_values(res_table)
    prob_INS.fillna(0.00, inplace=True)
    
    # Table of probabilities normalized by column:
    prob_INS_col_normalized = normalize_columns(prob_INS)

    # Get how many times each event is in the VCF table
    name_distribution = merged_df['name'].value_counts()

    # Transform it to a df
    name_distribution_df = pd.DataFrame({'name': name_distribution.index, 'number': name_distribution.values})

    # Calculate probabilities
    total = name_distribution_df['number'].sum()
    name_distribution_df['probability'] = name_distribution_df['number'] / total

    # # Generate the df with random numbers based on the probabilities
    sampled_names = np.random.choice(name_distribution_df['name'], size=num_events, p=name_distribution_df['probability'])
    df_insertions = pd.DataFrame({'name': sampled_names})

    # Add the info for the deletions
    df3 = add_columns(merged_df, prob_INS_col_normalized, df_insertions)

    # Add a new column 'Event_Type' with the value 'insertion'
    df3['Event_Type'] = 'Deletion'

    # Reorder columns to make 'Event_Type' the 4th column
    new_order = df3.columns.tolist()
    new_order.remove('Event_Type')
    new_order.insert(3, 'Event_Type')

    # Update order and name of columns
    df3 = df3[new_order]
    new_order = ['#ref', 'beg', 'end', 'Event_Type', 'name', 'Length']
    df_copy = df3[new_order].copy()
    df_copy.rename(columns={'Length': 'Total_Length'}, inplace=True)

    # Save final output
    df_copy.to_csv('Deletions_table.tsv', sep='\t', index=False)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Process deletions from VCF to BED format.')
    parser.add_argument('vcf_path_all', type=str, help='Path to the VCF file containing all deletions.')
    parser.add_argument('vcf_path_class', type=str, help='Path to the VCF file containing classified deletions.')
    parser.add_argument('chromosome_length', type=str, help='Path to the chromosome length file.')
    parser.add_argument('num_events', type=int, help='Number of events to sample (mandatory).')
    parser.add_argument('--bin_size', type=int, default=1000000, help='Size of genomic bins (default: 1000000).')

    args = parser.parse_args()
    main(args.vcf_path_all, args.vcf_path_class, args.bin_size, args.chromosome_length, args.num_events)
