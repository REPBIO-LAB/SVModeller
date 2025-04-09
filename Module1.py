# SVModeller - Module 1 

# Obtain and process data from Variant Calling File (VCF)

#Input:
# - VCF with insertion data
# - Chromosomes length (chr_length.txt)
# OPTIONAL: Window size for genome segmentation, by default 1 Mega base (integer number)

# Output: 
# - Insertion features (.tsv)
# - Genome-wide distribution (.tsv)
# - Probabilities (.tsv)
# - List of VNTR motifs (.txt)
# - List of SVA VNTR motifs (.txt)

import formats
import pandas as pd
import statistics
import gRanges
import argparse

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
            'Type_SV': variant.info.get('ITYPE_N', 'NA'),
            'Family': variant.info.get('FAM_N', 'NA'),
            'Conformation': variant.info.get('CONFORMATION', 'NA'),
            'Conformation_Ext': variant.info.get('CONFORMATION_EXT', 'NA'),
            'Start_position': variant.pos,
            'End_position': variant.pos + len(variant.alt),
            'Length': len(variant.alt),
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
            'Complete_Sequence': variant.alt
        })
    
    # Create a data frame from the extracted data
    df = pd.DataFrame(data, columns=[
        'Chromosome', 'Start_position', 'End_position', 'Type_SV', 'Family', 'Conformation', 'Conformation_Ext', 'Canonical', 
        'Length', 'PolyA_Length', 'Strand', 'TSD_Length', 'TD_5_Num', 'TD_5', 'TD_3_Num', 'TD_3', 'SVA_Hexamer', 'SVA_VNTR_Length', 
        'TD_orphan_Length', 'VNTR_Num_Motifs', 'VNTR_Motifs', 'SVA_VNTR_Coordinates', 'Complete_Sequence'
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

def SVA_VNTR_Motif(df):
    # Function to extract motif based on coordinates
    def extract_motif(row):
        coordinates = row['SVA_VNTR_Coordinates']
        sequence = row['Complete_Sequence']
        
        if pd.notna(coordinates) and '-' in coordinates:
            start, end = map(int, coordinates.split('-'))
            # Extract the segment from the sequence using the coordinates
            return sequence[start:end]
        return None  # Return None if there are no coordinates or the format is incorrect
    
    # Apply the extract_motif function to each row
    df['SVA_VNTR_Motif'] = df.apply(extract_motif, axis=1)
    
    # Drop the 'SVA_VNTR_Coordinates' and 'Complete_Sequence' columns
    df = df.drop(columns=['SVA_VNTR_Coordinates', 'Complete_Sequence'])
    
    return df

def extract_SVA_VNTR_Motifs(df):
    filename="SVA_VNTR_Motifs.txt"
    with open(filename, "w") as f:
        # Iterate through each row in the dataframe
        for motif in df['SVA_VNTR_Motif']:
            # Only write non-null motifs to the file
            if pd.notna(motif):
                f.write(str(motif) + "\n")
    
    # Drop the 'SVA_VNTR_Motif' column from the DataFrame
    df = df.drop(columns=['SVA_VNTR_Motif'])
    
    return df

def process_bed_table(result_df):
    # Change the names to the correct ones:
    result_BED_table = result_df.rename(columns={'Type_SV': 'name', 'Chromosome':'#ref', 'Start_position': 'beg', 'End_position':'end', 'Family':'SubType'})

    # Delete rows where 'name' is 'NA'
    result_BED_table = result_BED_table[result_BED_table['name'] != 'NA']
    
    # Delete rows where 'Canonical' is 'NOT_CANONICAL', except for those with 'name' as 'orphan'
    result_BED_table = result_BED_table[~((result_BED_table['Canonical'] == 'NOT_CANONICAL') & (result_BED_table['name'] != 'orphan'))]

    # Drop the 'Canonical' column after filtering
    result_BED_table.drop('Canonical', axis=1, inplace=True)

    # Delete rows where 'name' is 'DUP_INTERSPERSED'
    result_BED_table = result_BED_table[result_BED_table['name'] != 'DUP_INTERSPERSED']

    # Delete rows where 'name' is 'COMPLEX_DUP'
    result_BED_table = result_BED_table[result_BED_table['name'] != 'COMPLEX_DUP']
    
    # Extract the VNTR motifs from SVAs
    result_BED_table = SVA_VNTR_Motif(result_BED_table)
    result_BED_table = extract_SVA_VNTR_Motifs(result_BED_table)
    
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
    
    return result_BED_table

def create_dict(df):
    # Create the 'Event' column: Combine the column name and 'Conformation'
    df['Event'] = df.apply(lambda row: f"{row['name']}__{row['Conformation']}", axis=1)
    df['Event'] = df['Event'].replace(
            {'VNTR__NA': 'VNTR','DUP__NA': 'DUP', 'INV_DUP__NA': 'INV_DUP', 'NUMT__NA': 'NUMT', 'orphan__NA': 'orphan'}
        )
    # Create an empty dictionary to store results
    event_dict = {}

    # Iterate over each row in the DataFrame
    for index, row in df.iterrows():
        event = row['Event']  # Get the event value for the current row
        
        # Iterate over the columns from 'Length' (column 7) to the end
        for col in df.columns[6:]:  # Column 7 is index 7 (starting from 'Strand')
            value = row[col]  # Get the value in the current column
            if value != 'NA' and value != 'NA' and value is not None:  # Check if the value is not 'NA' or None
                # Construct the dictionary key
                key = f"{event}__{col}"
                
                # Add the value to the dictionary (create a list if key does not exist)
                if key not in event_dict:
                    event_dict[key] = []
                event_dict[key].append(value)

    # Remove keys that end with 'Event'
    keys_to_remove = [key for key in event_dict if key.endswith('Event')]
    
    # Delete the keys from the dictionary
    for key in keys_to_remove:
        del event_dict[key]

    return event_dict

def filter_sd(dict, key_list):
    ''' 
    Function to remove those values that exceed twice the standard deviation
    '''
    for key in key_list:
        values = dict.get(key, [])
        std_dev = statistics.stdev(values)
        dict[key] = [value for value in values if value <= 2 * std_dev]  
    return dict 

def extractvntr_motifs(dict_mutations):
    ''' 
    Extracts the motifs from the 'VNTR_VNTR_Motifs' key in the dictionary, 
    splits them into individual motifs, stores them in a list, 
    and then removes the 'VNTR_VNTR_Motifs' key from the dictionary.
    It also returns a DataFrame with the unique motifs and their proportions.
    
    Input:
    - dict_mutations (dict): The dictionary containing the VNTR data.
    
    Output:
    - dict_mutations (dict): The modified dictionary without the 'VNTR_VNTR_Motifs' key.
    - motifs_df (DataFrame): A DataFrame with unique motifs and their proportions.
    '''
    
    # Initialize the list for separated motifs
    separated_motifs = []

    # Check if 'VNTR_VNTR_Motifs' key exists in the dictionary
    if 'VNTR__VNTR_Motifs' in dict_mutations:
        # Get the list of values associated with the 'VNTR_VNTR_Motifs' key
        motifs = dict_mutations['VNTR__VNTR_Motifs']
        
        # Iterate through each motif and split by comma, extending the separated_motifs list
        for motif in motifs:
            separated_motifs.extend(motif.split(','))

        # Remove the 'VNTR_VNTR_Motifs' key from the dictionary
        del dict_mutations['VNTR__VNTR_Motifs']

    # Count occurrences of each motif in separated_motifs
    motif_counts = pd.Series(separated_motifs).value_counts()

    # Calculate proportions
    motif_proportions = motif_counts / motif_counts.sum()

    # Create a DataFrame with motifs and their proportions
    motifs_df = pd.DataFrame({
        'Motif': motif_proportions.index,
        'Proportion': motif_proportions.values
    })

    # Return the modified dictionary and the DataFrame
    return dict_mutations, motifs_df

def process_dictionary(dict_mutations):
    ''' 
    Process the dictionary by applying the filter_sd function for specific keys 
    and extracting VNTR motifs, then saving the motifs as a TSV file.
    
    Arguments:
    - dict_mutations (dict): The dictionary to process.
    
    Returns:
    - dict_mutations (dict): The modified dictionary after applying both functions.
    - motifs_df (DataFrame): A DataFrame containing the motifs and their proportions.
    '''
    # Step 1: Apply the filter_sd function to filter values based on standard deviation
    keys_filter = ['INV_DUP__Length', 'DUP__Length', 'VNTR__Length']
    dict_mutations = filter_sd(dict_mutations.copy(), keys_filter)
    
    # Step 2: Extract VNTR motifs and remove the 'VNTR_VNTR_Motifs' key
    dict_mutations, motifs_df = extractvntr_motifs(dict_mutations)
    
    # Step 3: Save the motifs DataFrame to a TSV file
    motifs_df.to_csv('Separated_Motifs.tsv', sep='\t', index=False)  # Save as TSV (tab-separated)
    
    # Return the modified dictionary and the DataFrame containing the motifs and their proportions
    return dict_mutations, motifs_df

def insertion_features_df(input_dict):
    # Define the columns as per the requirement
    columns = [
        'Event', 'Length', 'Strand', 'TSD_Length', 'TD_5', 'TD_3', 'SVA_Hexamer',
        'SVA_VNTR_Length', 'TD_orphan_Length', 'VNTR_Num_Motifs',
        'PolyA_Length_1', 'PolyA_Length_2', 'FOR', 'TRUN', 'REV', 'DEL', 'DUP'
    ]

    # Define the possible events that correspond to rows in the DataFrame
    event_column = [
        'SVA__Hexamer+Alu-like+VNTR+SINE-R+POLYA', 'Alu__FOR+POLYA',
        'Alu__TRUN+FOR+POLYA', 'SVA__Alu-like+VNTR+SINE-R+POLYA',
        'SVA__MAST2+VNTR+SINE-R+POLYA', 'L1__TRUN+FOR+POLYA+TD+POLYA',
        'L1__FOR+POLYA', 'SVA__VNTR+SINE-R+POLYA', 'L1__TRUN+FOR+POLYA',
        'SVA__VNTR+SINE-R+POLYA+TD+POLYA', 'SVA__SINE-R+POLYA',
        'SVA__TD+MAST2+VNTR+SINE-R+POLYA', 'orphan',
        'L1__TRUN+REV+DEL+FOR+POLYA', 'SVA__TD+Hexamer+Alu-like+VNTR+SINE-R+POLYA',
        'SVA__Hexamer+Alu-like+VNTR+SINE-R+POLYA+TD+POLYA',
        'L1__TD+FOR+POLYA', 'L1__TRUN+REV+DUP+FOR+POLYA',
        'SVA__Alu-like+VNTR+SINE-R+POLYA+TD+POLYA', 'L1__FOR+POLYA+TD+POLYA',
        'L1__TRUN+REV+DUP+FOR+POLYA+TD+POLYA',
        'L1__TRUN+REV+DEL+FOR+POLYA+TD+POLYA',
        'L1__TRUN+REV+BLUNT+FOR+POLYA', 'SVA__SINE-R+POLYA+TD+POLYA',
        'SVA__MAST2+VNTR+SINE-R+POLYA+TD+POLYA', 'L1__REV+DEL+FOR+POLYA',
        'L1__TRUN+REV+BLUNT+FOR+POLYA+TD+POLYA', 'VNTR', 'DUP',
        'INV_DUP', 'NUMT'
    ]

    # Step 1: Process the dictionary values to convert floats to integers
    for key, value in input_dict.items():
        # Convert floats to integers, leaving other types unchanged
        input_dict[key] = [
            (int(x) if isinstance(x, float) else x) for x in value
        ]
    
    # Create an empty DataFrame with the specified columns
    df = pd.DataFrame(columns=columns)
    
    # Step 2: Process the dictionary and add data to the DataFrame
    for key, values in input_dict.items():
        # Split the key into event name and column name (e.g., 'SVA_Hexamer+Alu-like+VNTR+SINE-R+POLYA_Length')
        event_name, column_name = key.rsplit('__', 1)

        # Check if the event name is valid (exists in the event_column list)
        if event_name in event_column and column_name in columns:
            # If event_name is not already in the DataFrame, create a new row with empty values
            if event_name not in df['Event'].values:
                # Create an empty row and add the event name
                empty_row = {col: '' for col in columns}
                empty_row['Event'] = event_name
                df = pd.concat([df, pd.DataFrame([empty_row])], ignore_index=True)
            
            # Get the index of the row that corresponds to the current event
            event_index = df[df['Event'] == event_name].index[0]
            
            # Assign the values to the correct cell, join by commas in case of multiple values
            # Avoid overwriting cells if they already have values
            existing_value = df.at[event_index, column_name]
            if existing_value:
                new_value = ','.join(map(str, values)) if values else ''
                # Check for duplicates and join without repeating values
                existing_value = set(existing_value.split(','))
                new_value = set(new_value.split(','))
                final_value = list(existing_value | new_value)
                df.at[event_index, column_name] = ','.join(final_value)
            else:
                df.at[event_index, column_name] = ','.join(map(str, values)) if values else ''
    
    # Save the Insertion Features df to a .tsv file
    df.to_csv('Insertion_Features.tsv', sep='\t', index=False)

def probabilities_df(table):
    # Get how many times each event is in the VCF table
    name_distribution = table['Event'].value_counts()
    
    # Transform it to a df
    name_distribution_df = pd.DataFrame({'Event': name_distribution.index, 'number': name_distribution.values})

    # Calculate probabilities
    total = name_distribution_df['number'].sum()
    name_distribution_df['Probability'] = name_distribution_df['number'] / total
    name_distribution_df = name_distribution_df.drop(columns=['number'])

    name_distribution_df.to_csv('Probabilities.tsv', sep='\t', index=False)

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
    event_column = [
        'SVA__Hexamer+Alu-like+VNTR+SINE-R+POLYA', 'Alu__FOR+POLYA',
        'Alu__TRUN+FOR+POLYA', 'SVA__Alu-like+VNTR+SINE-R+POLYA',
        'SVA__MAST2+VNTR+SINE-R+POLYA', 'L1__TRUN+FOR+POLYA+TD+POLYA',
        'L1__FOR+POLYA', 'SVA__VNTR+SINE-R+POLYA', 'L1__TRUN+FOR+POLYA',
        'SVA__VNTR+SINE-R+POLYA+TD+POLYA', 'SVA__SINE-R+POLYA',
        'SVA__TD+MAST2+VNTR+SINE-R+POLYA', 'orphan',
        'L1__TRUN+REV+DEL+FOR+POLYA', 'SVA__TD+Hexamer+Alu-like+VNTR+SINE-R+POLYA',
        'SVA__Hexamer+Alu-like+VNTR+SINE-R+POLYA+TD+POLYA',
        'L1__TD+FOR+POLYA', 'L1__TRUN+REV+DUP+FOR+POLYA',
        'SVA__Alu-like+VNTR+SINE-R+POLYA+TD+POLYA', 'L1__FOR+POLYA+TD+POLYA',
        'L1__TRUN+REV+DUP+FOR+POLYA+TD+POLYA',
        'L1__TRUN+REV+DEL+FOR+POLYA+TD+POLYA',
        'L1__TRUN+REV+BLUNT+FOR+POLYA', 'SVA__SINE-R+POLYA+TD+POLYA',
        'SVA__MAST2+VNTR+SINE-R+POLYA+TD+POLYA', 'L1__REV+DEL+FOR+POLYA',
        'L1__TRUN+REV+BLUNT+FOR+POLYA+TD+POLYA', 'VNTR', 'DUP',
        'INV_DUP', 'NUMT'
    ]

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


def genome_wide_distribution(chromosome_length, bin_size, table):
    # Dictionary containing references as keys and their lengths as values:
    chr_length = formats.chrom_lengths_index(chromosome_length)

    # Bin size
    binSize = bin_size  # this is 1MB
    
    # Generate the list of chromosomes  (assuming 'chr1' to 'chr22')
    chromosomes = [f'chr{i}' for i in range(1, 23)]
    
    # Create genomic bins based on chromosome lengths and bin size
    bins = gRanges.makeGenomicBins(chr_length, binSize, chromosomes)[::-1]
    
    # Create the table of insertions classified in windows
    res_table = mut_bins(bins, table)
    final_table = normalize_columns(res_table)

    # Save the result to a TSV file
    final_table.to_csv('Genome_Wide_Distribution.tsv', sep='\t', index=False)


def main(file_path, chromosome_length, bin_size):
    print(f'VCF file with insertions: {file_path}')
    print(f'File with chromosomes length: {chromosome_length}')
    print(f'Size of genomic bins (default: 1000000): {bin_size}')

    table = read_vcf_file_BED(file_path)
    processed_table = process_bed_table(table)
    dict1 = create_dict(processed_table)
    processed_dict, motifs_df = process_dictionary(dict1)
    insertion_features_df(processed_dict)
    genome_wide_distribution(chromosome_length,bin_size,processed_table)
    probabilities_df(processed_table)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Obtain and process data from Variant Calling File (VCF)')
    parser.add_argument('file_path', type=str, help='Path to the VCF file containing insertion data.')
    parser.add_argument('chromosome_length', type=str, help='Path to the chromosome length file.')
    parser.add_argument('--bin_size', type=int, default=1000000, help='Size of genomic bins (default: 1000000).')
    args = parser.parse_args()
    main(args.file_path, args.chromosome_length, args.bin_size)
