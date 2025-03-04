# SVMoldeler - Module 3

# Process Deletions

# INPUT:
# - VCF file with deletion data
# - Chromosome length file
#- Number of events
# - Size of genomic bins (default: 1000000)

# OTPUT:
# - Deletion events table (.tsv)

import argparse
import formats
import gRanges
import numpy as np
import pandas as pd

def TD_filter(df):
    ''' 
    Input & output: df
    1- checks if TD_5_Num is not NA and has a value higher than 1
    2- if so, deletes the value of TD_5_Num and TD_5
    3- does the same for TD_3_Num
    4- deletes columns TD_5_Num and TD_3_Num at the end
    '''
    # Iterate over each row in the DataFrame
    for index, row in df.iterrows():
        # Check and modify TD_5_Num and TD_5 if condition is met
        if row['TD_5_Num'] != 'NA' and int(row['TD_5_Num']) > 1:
            df.at[index, 'TD_5'] = 'NA'
            df.at[index, 'TD_5_Num'] = 'NA'
        
        # Check and modify TD_3_Num and TD_3 if condition is met
        if row['TD_3_Num'] != 'NA' and int(row['TD_3_Num']) > 1:
            df.at[index, 'TD_3'] = 'NA'
            df.at[index, 'TD_3_Num'] = 'NA'
    
    # After processing, drop the TD_5_Num and TD_3_Num columns
    df = df.drop(columns=['TD_5_Num', 'TD_3_Num'])
    
    return df

def read_vcf_file_BED(file_path):
    ''' 
    Function to:
    - Initialize a VCF object
    - Use the VCF class to read the path to the VCF file
    - Create an empty list to store the data
    - For each line, get the information and place it in a new DataFrame
    - First, the values that are not common in all lines are defined to get their value or write NA
    - Add to the list all the information and add this list to the DataFrame 
    '''
    vcf = formats.VCF()  # Instantiate the VCF class
    vcf.read(file_path)  # Read the VCF file

    data = []  # List to store the extracted data
    
    for variant in vcf.variants:
        # Check for "NOT_CANONICAL" in the INFO field
        canonical_status = 'NOT_CANONICAL' if 'NOT_CANONICAL' in variant.info else ''

        # Extract information from the INFO field, and if not present, set it to NA
        data.append({
            'Type_SV': variant.info.get('DTYPE_N', 'NA'),
            'Family': variant.info.get('FAM_N', 'NA'),
            'Conformation': variant.info.get('CONFORMATION', 'NA'),
            'Conformation_Ext': variant.info.get('CONFORMATION_EXT', 'NA'),
            'Start_position': variant.pos,
            'End_position': variant.pos + len(variant.alt),
            'Length': variant.info.get('DEL_LEN', 'NA'),
            'Chromosome': variant.chrom,
            'PolyA_Length': variant.info.get('POLYA_LEN', 'NA'),
            'Strand': variant.info.get('STRAND', 'NA'),
            'TSD_Length': variant.info.get('TSD_LEN', 'NA'),
            'TD_5_Num': variant.info.get('5PRIME_NB_TD', 'NA'),
            'TD_5': variant.info.get('5PRIME_TD_LEN', 'NA'),
            'TD_3_Num': variant.info.get('3PRIME_NB_TD', 'NA'),
            'TD_3': variant.info.get('3PRIME_TD_LEN', 'NA'),
            'SVA_Hexamer': variant.info.get('HEXAMER_LEN', 'NA'),
            'TD_orphan_Length': variant.info.get('ORPHAN_TD_LEN', 'NA'),
            'VNTR_Num_Motifs': variant.info.get('NB_MOTIFS', 'NA'),
            'VNTR_Motifs': variant.info.get('MOTIFS', 'NA'),
            'Canonical': canonical_status,
            'SVA_VNTR_Length': variant.info.get('VNTR_LEN', 'NA'),
            'SVA_VNTR_Coordinates': variant.info.get('VNTR_COORD', 'NA'),
        })
    
    # Create a data frame from the extracted data
    df = pd.DataFrame(data, columns=[
        'Chromosome', 'Start_position', 'End_position', 'Type_SV', 'Family', 'Conformation', 'Conformation_Ext', 'Canonical', 
        'Length', 'PolyA_Length', 'Strand', 'TSD_Length', 'TD_5_Num', 'TD_5', 'TD_3_Num', 'TD_3', 'SVA_Hexamer', 'SVA_VNTR_Length', 
        'TD_orphan_Length', 'VNTR_Num_Motifs', 'VNTR_Motifs', 'SVA_VNTR_Coordinates'
    ])
    
    return df

def extract_conformation_data(row):
    conformation_ext = row['Conformation_Ext']
    
    # Initialize columns with 'NA'
    for_val = trun_val = rev_val = del_val = dup_val = 'NA'
    
    # Extract the data before the commas in the Conformation_Ext column
    if 'FOR' in conformation_ext:
        # Get the number before the comma in 'FOR'
        for_val = conformation_ext.split('FOR(')[1].split(',')[0] if 'FOR' in conformation_ext else 'NA'
    if 'TRUN' in conformation_ext:
        # Get the number before the comma in 'TRUN'
        trun_val = conformation_ext.split('TRUN(')[1].split(',')[0] if 'TRUN' in conformation_ext else 'NA'
    if 'REV' in conformation_ext:
        # Get the number before the comma in 'REV'
        rev_val = conformation_ext.split('REV(')[1].split(',')[0] if 'REV' in conformation_ext else 'NA'
    if 'DEL' in conformation_ext:
        # DEL doesn't have a comma, we just extract the value inside DEL()
        del_val = conformation_ext.split('DEL(')[1].split(')')[0] if 'DEL' in conformation_ext else 'NA'
    if 'DUP' in conformation_ext:
        # DUP doesn't have a comma, we just extract the value inside DEL()
        dup_val = conformation_ext.split('DUP(')[1].split(')')[0] if 'DUP' in conformation_ext else 'NA'   

    return pd.Series([for_val, trun_val, rev_val, del_val, dup_val])

def process_bed_table(result_df):
    # Change the names to the correct ones:
    result_BED_table = result_df.rename(columns={'Type_SV': 'name', 'Chromosome':'#ref', 'Start_position': 'beg', 'End_position':'end', 'Family':'SubType'})

    # Update 'name' column where 'name' is 'NA' to 'simple'
    result_BED_table['name'] = result_BED_table['name'].replace('NA', 'simple')

    # Delete rows where 'Canonical' is 'NOT_CANONICAL', except for those with 'name' as 'orphan'
    result_BED_table = result_BED_table[~((result_BED_table['Canonical'] == 'NOT_CANONICAL') & (result_BED_table['name'] != 'orphan'))]

    # Drop the 'Canonical' column after filtering
    result_BED_table.drop('Canonical', axis=1, inplace=True)

    # Delete rows where 'name' is 'DUP_INTERSPERSED'
    result_BED_table = result_BED_table[result_BED_table['name'] != 'DUP_INTERSPERSED']

    # Delete rows where 'name' is 'COMPLEX_DUP'
    result_BED_table = result_BED_table[result_BED_table['name'] != 'COMPLEX_DUP']
    
    # Update 'name' for rows with 'solo' or 'partnered' based on 'SubType'
    result_BED_table['name'] = result_BED_table.apply(
        lambda row: row['SubType'] if row['name'] in ['solo', 'partnered'] else row['name'], axis=1
    )

    # Delete the 'SubType' column
    result_BED_table = result_BED_table.drop('SubType', axis=1)

    # Applying the filter to get only the ones that have 1 transduction
    result_BED_table = TD_filter(result_BED_table)

    # From the PolyA_Length, get the 1st and 2nd values in separate columns
    result_BED_table['PolyA_Length_1'] = result_BED_table['PolyA_Length'].apply(lambda x: x.split(',')[0] if x != 'NA' else x)
    result_BED_table['PolyA_Length_2'] = result_BED_table['PolyA_Length'].apply(lambda x: x.split(',')[1] if x != 'NA' and len(x.split(',')) > 1 else 'NA')

    # Delete the original PolyA_Length column
    result_BED_table.drop('PolyA_Length', axis=1, inplace=True)

    # Convert 'PolyA_Length' columns to numeric, converting 'NA' to NaN
    result_BED_table['PolyA_Length_1'] = pd.to_numeric(result_BED_table['PolyA_Length_1'], errors='coerce')
    result_BED_table['PolyA_Length_2'] = pd.to_numeric(result_BED_table['PolyA_Length_2'], errors='coerce')

    # Transform TSD_Length & SVA_Hexamer columns to numeric
    result_BED_table['TSD_Length'] = pd.to_numeric(result_BED_table['TSD_Length'], errors='coerce')
    result_BED_table['SVA_Hexamer'] = pd.to_numeric(result_BED_table['SVA_Hexamer'], errors='coerce')
    
    # Replace all NaN values with 'NA'
    result_BED_table.fillna('NA', inplace=True)
    
    # Apply the extraction function to create the new columns
    result_BED_table[['FOR', 'TRUN', 'REV', 'DEL', 'DUP']] = result_BED_table.apply(extract_conformation_data, axis=1)
    
    # Create the 'Event' column: Combine the column name and 'Conformation'
    result_BED_table['Event'] = result_BED_table.apply(lambda row: f"{row['name']}__{row['Conformation']}", axis=1)
    result_BED_table['Event'] = result_BED_table['Event'].replace(
            {'VNTR__NA': 'VNTR','DUP__NA': 'DUP', 'INV_DUP__NA': 'INV_DUP', 'NUMT__NA': 'NUMT', 'orphan__NA': 'orphan', 'simple__NA': 'simple'}
        )
    
    # Delete the 'SubType' column
    result_BED_table = result_BED_table.drop('Conformation_Ext', axis=1)
    return result_BED_table

def mut_bins(bins, table):
    """
    Classify mutations in each window of the bindatabase.

    Parameters:
    bins (list): List of windows in the format [(chromosome, start, end), ...]
    table (pandas.DataFrame): Table with mutation data

    Returns:
    pandas.DataFrame: Classified mutations in each window
    """
    # Create a dictionary to map event types to column names
    event_column = ['simple', 'Alu__TRUN+FOR+POLYA', 'orphan', 'L1__TRUN+FOR+POLYA',
       'SVA__TD+MAST2+VNTR+SINE-R+POLYA', 'SVA__MAST2+VNTR+SINE-R+POLYA',
       'L1__FOR+POLYA', 'L1__TD+FOR+POLYA',
       'SVA__Hexamer+Alu-like+VNTR+SINE-R+POLYA',
       'L1__TRUN+REV+DEL+FOR+POLYA', 'Alu__FOR+POLYA',
       'SVA__VNTR+SINE-R+POLYA', 'L1__TRUN+REV+DEL+FOR+POLYA+TD+POLYA',
       'SVA__Alu-like+VNTR+SINE-R+POLYA', 'SVA__SINE-R+POLYA',
       'SVA__Hexamer+Alu-like+VNTR+SINE-R+POLYA+TD+POLYA',
       'L1__TRUN+FOR+POLYA+TD+POLYA',
       'SVA__MAST2+VNTR+SINE-R+POLYA+TD+POLYA',
       'L1__TRUN+REV+DUP+FOR+POLYA', 'VNTR']

    # Initialize the output DataFrame with columns: 'window', 'beg', 'end', and the event types
    df = pd.DataFrame(columns=['window', 'beg', 'end'] + event_column)

    # Iterate over each window in bins
    for window in bins:
        chrom, start, end = window
        # Filter mutations that fall within the current window
        window_df = table[(table['#ref'] == chrom) &
                                     (table['beg'] >= start) &
                                     (table['end'] <= end)]

        # Count the occurrences of each event type in the window
        counts = window_df['Event'].value_counts().to_dict()

        # Create a row for the current window
        row = [chrom, start, end]

        # For each event type, append the count (or 0 if not present) to the row
        for event in event_column:
            row.append(counts.get(event, 0))

        # Append the row to the output DataFrame
        df.loc[len(df)] = row

    return df

def classify_mutations_in_bins(chromosome_length, bin_size, merged_df):
    """
    Classify mutations into genomic bins based on chromosome length and bin size.

    Parameters:
    chromosome_length (dict): Dictionary with chromosome lengths.
    bin_size (int): The size of each genomic bin.
    merged_df: DataFrame containing mutation data.

    Returns:
    pandas.DataFrame: Classified mutations for each genomic bin.
    """
    # Get chromosome lengths
    chr_length = formats.chrom_lengths_index(chromosome_length)

    # Generate bins for the genome, including 'chrX'
    bins = gRanges.makeGenomicBins(chr_length, bin_size, [f'chr{i}' for i in range(1, 23)] + ['chrX'])[::-1]

    # Classify mutations in each window of the genomic bins
    res_table = mut_bins(bins, merged_df)

    return res_table

def normalize_columns(df):
    '''
    Normalize columns in the dataframe starting from the 3rd column to the end.
    Each column will be normalized by dividing by the sum of the column values.
    '''
    # Select columns from the 3rd column onward
    columns_to_normalize = df.columns[3:]
    
    for column in columns_to_normalize:
        # Get the sum of the column
        total_sum = df[column].sum()
    
        # If the total sum is not zero, normalize the column
        if total_sum != 0:
            df[column] = df[column] / total_sum
    
    return df

def probabilities_df(table):
    # Get how many times each event is in the VCF table
    name_distribution = table['Event'].value_counts()
    
    # Transform it to a df
    name_distribution_df = pd.DataFrame({'Event': name_distribution.index, 'number': name_distribution.values})

    # Calculate probabilities
    total = name_distribution_df['number'].sum()
    name_distribution_df['Probability'] = name_distribution_df['number'] / total
    
    # Drop the 'number' column
    name_distribution_df = name_distribution_df.drop(columns=['number'])
    
    return name_distribution_df

def add_columns(df1, df2, df3):
    ''' 
    Function to add start and end columns to the df based on probabilities
    '''
    new_df = pd.DataFrame(columns=['#ref', 'beg', 'end', 'Length'])
    for _, row in df3.iterrows():
        name = row['Event']
        prob_df = df2[df2[name] > 0]
        if prob_df.empty:
            random_row = df1[df1['Event'] == name].sample(n=1)
        else:
            name_df = df1[df1['Event'] == name]
            if name_df.empty:
                continue
            weights = prob_df[name].values
            if weights.sum() == 0:
                weights = np.ones(len(weights)) / len(weights)
            weights = np.repeat(weights, len(name_df) // len(weights) + 1)[:len(name_df)]
            random_row = name_df.sample(n=1, weights=weights)
        new_df = new_df._append(random_row[['#ref', 'beg', 'end', 'Length']], ignore_index=True)
    df3[['#ref', 'beg', 'end', 'Length']] = new_df
    return df3

def generate_deletion_events(probabilities_table, num_events, deletions_table, genome_wide_distribution):
    # Generate name of the events based on their proportions
    sampled_names = np.random.choice(probabilities_table['Event'], size=num_events, p=probabilities_table['Probability'])
    table_events = pd.DataFrame({'Event': sampled_names})

    # Add the info for the deletions
    df3 = add_columns(deletions_table, genome_wide_distribution, table_events)
    
    # Add a new column 'Event_Type' with the value 'insertion'
    df3['Event_Type'] = 'Deletion'

    # Define the new order of columns first
    new_order = ['#ref', 'beg', 'end', 'Event_Type', 'Event', 'Length']
    
    # Update order and name of columns
    df3 = df3[new_order]
    
    # Copy the DataFrame and rename 'Length' column to 'Total_Length'
    df_copy = df3.copy()
    df_copy.rename(columns={'Length': 'Total_Length'}, inplace=True)

    return df_copy

def main(vcf_path, path_chromosome_length, num_events, bin_size,):
    # Print the paths of the input files
    print(f'VCF file with deletion data: {vcf_path}')
    print(f'Chromosome length file: {path_chromosome_length}')
    print(f'Number of events: {num_events}')
    print(f'Size of genomic bins (default: 1000000).: {bin_size}')

    # Get data from VCF file
    table = read_vcf_file_BED(vcf_path)
    
    # Process the data
    processed_table = process_bed_table(table)
    
    # Obtain genome-wide distribution of the events & normalize it
    genome_wide_distribution = classify_mutations_in_bins(path_chromosome_length,bin_size,processed_table)
    normalize_columns(genome_wide_distribution)
    
    # Obtain probability of each event
    probabilities = probabilities_df(processed_table)
    
    # Generate the deletion events & save the output
    deletion_events = generate_deletion_events(probabilities,num_events,processed_table,genome_wide_distribution)
    deletion_events.to_csv('Deletions_table.tsv', sep='\t', index=False)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Process deletions from VCF to BED format.')
    parser.add_argument('vcf_path', type=str, help='Path to the VCF file containing deletion data.')
    parser.add_argument('path_chromosome_length', type=str, help='Path to the chromosome length file.')
    parser.add_argument('num_events', type=int, help='Number of events to sample (mandatory).')
    parser.add_argument('--bin_size', type=int, default=1000000, help='Size of genomic bins (default: 1000000).')

    args = parser.parse_args()
    main(args.vcf_path, args.path_chromosome_length, args.num_events, args.bin_size,)
