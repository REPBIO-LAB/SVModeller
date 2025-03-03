# SVModeler - Module 2

# Model probabilities, insertion features, and genome wide distribution to generate new events

import argparse
import numpy as np
import pandas as pd
import distfit
from distfit import distfit
import random
import statistics
import pysam

def consensus_seqs(file_path):
    ''' 
    Function to read from a fasta file the consensus sequences

    Input: path of the fasta file
    Output: dictionary containing, for example, key: Alu_Seq - Value: the consensus sequence
    '''
    sequences = {"Alu_Seq": "", "L1_Seq": "", "SVA_Alu-like_Seq": "", "SVA_SINE-R_Seq": "",
                "SVA_MAST2_Seq": "", "NUMT_Seq": ""}
    
    with open(file_path, "r") as file:
        lines = file.readlines()
        
    current_sequence = "" 
    previous_header = ""  # Initialize previous_header here to avoid the UnboundLocalError

    for i in range(len(lines)):
        line = lines[i].strip()  # Remove any leading/trailing spaces or newlines
        
        if line.startswith(">"): 
            # If we were already collecting a sequence, assign it to the corresponding key
            if current_sequence:
                if "consensus|Alu" in previous_header:
                    sequences["Alu_Seq"] = current_sequence
                elif "consensus|L1" in previous_header:
                    sequences["L1_Seq"] = current_sequence
                elif "consensus|SVA|SVA_F|Alu-like" in previous_header:
                    sequences["SVA_Alu-like_Seq"] = current_sequence
                elif "consensus|SVA|SVA_F|SINE-R" in previous_header:
                    sequences["SVA_SINE-R_Seq"] = current_sequence
                elif "consensus|SVA|exon1|MAST2" in previous_header:
                    sequences["SVA_MAST2_Seq"] = current_sequence
                elif "consensus|NC_012920.1" in previous_header:
                    sequences["NUMT_Seq"] = current_sequence
            # Reset current sequence and store the new header
            current_sequence = ""
            previous_header = line  # Keep track of the header to check for matching keys
        else:
            # Append the current line (sequence part) to the current sequence
            current_sequence += line.strip()  # Remove extra spaces or newlines
    
    # Don't forget to handle the last sequence
    if current_sequence:
        if "consensus|Alu" in previous_header:
            sequences["Alu_Seq"] = current_sequence
        elif "consensus|L1" in previous_header:
            sequences["L1_Seq"] = current_sequence
        elif "consensus|SVA|SVA_F|Alu-like" in previous_header:
            sequences["SVA_Alu-like_Seq"] = current_sequence
        elif "consensus|SVA|SVA_F|SINE-R" in previous_header:
            sequences["SVA_SINE-R_Seq"] = current_sequence
        elif "consensus|SVA|exon1|MAST2" in previous_header:
            sequences["SVA_MAST2_Seq"] = current_sequence
        elif "consensus|NC_012920.1" in previous_header:
            sequences["NUMT_Seq"] = current_sequence

    return sequences

def probabilities_total_number(probabilities_numbers_df,num_events):
    table = pd.read_csv(probabilities_numbers_df, sep='\t')
    print(f"Columns in the table: {table.columns}")
    if 'Probability' in table.columns:
        sampled_names = np.random.choice(table['Event'], size=num_events, p=table['Probability'])
        table_events = pd.DataFrame({'name': sampled_names})
    elif 'Number' in table.columns:
        rows = []
        for i, row in table.iterrows():
            event = row['Event']
            count = row['Number']
            rows.extend([event] * int(count))
        table_events = pd.DataFrame({'name': rows})
    else:
        raise ValueError("The second column must be either 'Probability' or 'Number'.")
    
    return table_events

def generate_dict_from_table(data):
    result_dict = {}
    
    for index, row in data.iterrows():
        event = row['Event']
        
        for column in data.columns[1:]:  # Skip the 'Event' column
            value = row[column]
            
            if pd.notna(value):  # Check if value is not NaN
                key = f"{event}__{column}"
                result_dict[key] = [int(x) for x in value.split(',') if x.isdigit()]
    
    # Remove keys that end with 'Strand' and 'Length'
    result_dict = {key: value for key, value in result_dict.items() if not key.endswith('Strand')}

    return result_dict

def generate_random_numbers(dictionary, num_samples):
    random_numbers_dict = {}
    for key, value in dictionary.items():
        dist = distfit()
        dist.fit_transform(value)
        random_numbers = dist.generate(num_samples)
        random_numbers = random_numbers.astype(int)
        random_numbers_dict[key] = random_numbers
    return random_numbers_dict

def remove_negative_values(values):
    ''' 
    Function to remove negative values from the list
    '''
    return [value for value in values if value > 0]

def filter_sd(dict,key_list):
    ''' 
    Function to remove those values that exceed twice the standard deviation
    '''
    for key in key_list:
        values = dict[key]
        std_dev = statistics.stdev(values)
        values = [value for value in values if value <= 2 * std_dev]
        dict[key] = values

def filter_DEL(dict):
    '''
    Check if a key contains 'DEL' in its name. If it does, it removes values larger than 1000 from the corresponding list of values for that key.
    '''
    for key in dict:
        if 'DEL' in key:
            # Get the list of values for this key
            values = dict[key]
            # Filter out values greater than 1000
            filtered_values = [value for value in values if value <= 1000]
            # Update the dictionary with the filtered values
            dict[key] = filtered_values
    return dict

def filter_FOR(dict_insertion_features, dict_consensus):
    # Extract the Alu_Seq and L1_Seq from dict_consensus and get their lengths
    aluseq_seq = dict_consensus.get('Alu_Seq', '')
    l1seq_seq = dict_consensus.get('L1_Seq', '')
    
    aluseq_length = len(aluseq_seq) if aluseq_seq else 0
    l1seq_length = len(l1seq_seq) if l1seq_seq else 0
    
    # Iterate over the dictionary
    for key, values in dict_insertion_features.items():
        if 'FOR' in key and key.endswith('__FOR'):  # Check if 'FOR' is in key and ends with '__FOR'
            if 'L1' in key:
                # Remove values larger than L1_seq length for 'FOR' and 'L1' keys
                dict_insertion_features[key] = [value for value in values if value <= l1seq_length]
            elif 'Alu' in key:
                # Remove values larger than Alu_seq length for 'FOR' and 'Alu' keys
                dict_insertion_features[key] = [value for value in values if value <= aluseq_length]

    return dict_insertion_features

def distribution_random_numbers(dict_insertion_features,num_events,dict_consensus):
    array_dict = {}
    for key in dict_insertion_features:
        # create a numpy array for each key
        array = np.array(dict_insertion_features[key])
        array_dict[key] = array

    # From the dictionary, generate the distributions, and from these generate a dictionary with determined random numbers based on the distributions.
    number_random_events=num_events*3
    random_numbers_dict = generate_random_numbers(array_dict, number_random_events)

    # Filter to remove possible negative values
    for key, values in random_numbers_dict.items():
        random_numbers_dict[key] = remove_negative_values(values)

    # Filter to remove possible larger values that exceed twice the standard deviation
    key_list = ['DUP__Length', 'INV_DUP__Length', 'VNTR__Length', 'NUMT__Length']
    filter_sd(random_numbers_dict,key_list)
    
    # Filter largest DEL values
    filter_DEL(random_numbers_dict)
    
    # Filter FOR values larger than the consensus of L1 and Alu
    filter_FOR(random_numbers_dict, dict_consensus)

    return random_numbers_dict

def process_insertion_features_random_numbers(insertion_features_df,num_events,dict_consensus):
    # Open insertion features df
    table_insertion_features = pd.read_csv(insertion_features_df, sep='\t')
    # Create dictionary from the df
    dict_insertion_features = generate_dict_from_table(table_insertion_features)
    # Generate distributions of data (disfit) and dictionary of random numbers for every event and feature
    dict_random = distribution_random_numbers(dict_insertion_features,num_events,dict_consensus)
    return dict_random

def add_beg_end_columns(df_insertions, genome_wide_distribution_df):
    ''' 
    Function to add ref and beg columns to the df of insertions based on genome-wide distribution
    '''
    # Open insertion features df
    genome_wide_distribution = pd.read_csv(genome_wide_distribution_df, sep='\t')

    # Create new columns in the first DataFrame
    df_insertions['#ref'] = ''
    df_insertions['beg'] = ''

    # Function to select a random row based on probabilities
    def select_random_row(probabilities):
        return np.random.choice(probabilities.index, p=probabilities)

    # Iterate over each row in the first DataFrame
    for index, row in df_insertions.iterrows():
        event_name = row['name']
        probabilities = genome_wide_distribution[event_name]
        selected_row = select_random_row(probabilities)
        
        # Fill the values in the first DataFrame
        df_insertions.at[index, '#ref'] = genome_wide_distribution.at[selected_row, 'window']
        df_insertions.at[index, 'beg'] = np.random.randint(genome_wide_distribution.at[selected_row, 'beg'], genome_wide_distribution.at[selected_row, 'end'])

    return df_insertions

def add_elements_columns(dict, df_insertions):
   
    # Create new columns with default value 0
    columns = ['Length', 'PolyA_Length_1', 'PolyA_Length_2', 'TSD_Length', 'TD_5', 'TD_3', 'TD_orphan_Length', 
               'VNTR_Num_Motifs', 'SVA_Hexamer', 'SVA_VNTR_Length', 'FOR', 'TRUN', 'REV', 'DEL', 'DUP']
    
    for column in columns:
        df_insertions[column] = 0

    # Iterate over rows in df_insertions
    for index, row in df_insertions.iterrows():
        name = row['name']
        for column in columns:
            key = f"{name}__{column}"
            # Check if the key exists in the dictionary and has corresponding values
            if key in dict and dict[key]:
                # Select a random value from the list without modifying the dictionary
                df_insertions.loc[index, column] = random.choice(dict[key])
            else:
                # If the key is not found or the list is empty, assign 0
                df_insertions.loc[index, column] = 0

    # Fill any NaN values with 0
    df_insertions.fillna(0, inplace=True)

    return df_insertions

def add_SVA_info(df, dict_consensus):
    # Extract the sequence data from dict_consensus
    aluseq_seq = dict_consensus.get('SVA_Alu-like_Seq', '')
    siner_seq = dict_consensus.get('SVA_SINE-R_Seq', '')
    mast2_seq = dict_consensus.get('SVA_MAST2_Seq', '')

    # Get the lengths of the sequences
    aluseq_length = len(aluseq_seq) if aluseq_seq else 0
    siner_length = len(siner_seq) if siner_seq else 0
    mast2_length = len(mast2_seq) if mast2_seq else 0

    # Create the new columns in the dataframe
    df['SINE_R'] = 0
    df['MAST2'] = 0
    df['ALU_LIKE'] = 0
    
    # Iterate through each row and check the 'name' column
    for index, row in df.iterrows():
        if 'Alu-like' in row['name']:
            df.at[index, 'ALU_LIKE'] = aluseq_length
        if 'SINE-R' in row['name']:
            df.at[index, 'SINE_R'] = siner_length
        if 'MAST2' in row['name']:
            df.at[index, 'MAST2'] = mast2_length
    
    return df

# Function to calculate the Length for each row
def calculate_length(row, columns_to_sum):
    # Check if 'name' is one of the excluded values
    if row['name'] in ['NUMT', 'DUP', 'INV_DUP', 'VNTR']:
        return row['Length']  # Do not change the Length for these rows
    
    # Sum all the relevant columns, considering missing values as 0
    total_sum = 0
    for col in columns_to_sum:
        # Ensure the column value is numeric, convert if necessary
        value = pd.to_numeric(row.get(col, 0), errors='coerce')  # Convert to numeric, coercing errors to NaN
        total_sum += value if not pd.isna(value) else 0  # Add the value or 0 if NaN
    
    # Subtract the value in the 'DEL' column, if it exists
    del_value = pd.to_numeric(row.get('DEL', 0), errors='coerce')
    total_sum -= del_value if not pd.isna(del_value) else 0  # Subtract the value or 0 if NaN
    
    return total_sum

def update_trun(df_inertions, dict_consensus):
    # Extract the Alu_Seq and L1_Seq from dict_consensus and get their lengths
    aluseq_seq = dict_consensus.get('Alu_Seq', '')
    l1seq_seq = dict_consensus.get('L1_Seq', '')
    
    # Get the length of the DNA sequences (length of the string)
    aluseq_length = len(aluseq_seq) if aluseq_seq else 0
    l1seq_length = len(l1seq_seq) if l1seq_seq else 0

    # Iterate over the DataFrame rows
    for index, row in df_inertions.iterrows():
        name_value = str(row['name'])
        
        # Calculate the sum of 'FOR', 'DEL', and 'REV' values, treating NaN as 0
        total_subtraction = 0
        for col in ['FOR', 'DEL', 'REV']:
            value = row.get(col, 0)
            
            # Ensure the value is numeric (if not, treat it as 0)
            try:
                value = float(value)
            except ValueError:
                value = 0  # If it's a non-numeric value, treat it as 0
            
            total_subtraction += value
        
        # Check if both 'L1' or 'Alu' and 'TRUN' are in the 'name' column and update 'TRUN' accordingly
        if 'L1' in name_value and 'TRUN' in name_value:
            # Calculate TRUN for 'L1'
            df_inertions.at[index, 'TRUN'] = l1seq_length - total_subtraction
        elif 'Alu' in name_value and 'TRUN' in name_value:
            # Calculate TRUN for 'Alu'
            df_inertions.at[index, 'TRUN'] = aluseq_length - total_subtraction

    return df_inertions

# Function to update the DataFrame
def update_dataframe(df_insertions, dict_consensus):
    # Add a column 'Strand' with randomly assigned '+' or '-'
    df_insertions['Strand'] = np.random.choice(['+', '-'], size=len(df_insertions))
    
    # List of columns to sum
    columns_to_sum = [
        'PolyA_Length_1', 'PolyA_Length_2', 'TSD_Length', 'TD_5', 'TD_3', 'TD_orphan_Length', 
        'SVA_Hexamer', 'SVA_VNTR_Length', 'FOR', 'REV', 'DUP', 'SINE_R', 'MAST2', 'ALU_LIKE'
    ]
    
    # Apply the calculate_length function to each row and update the 'Length' column
    df_insertions['Length'] = df_insertions.apply(calculate_length, axis=1, columns_to_sum=columns_to_sum)
    
    # Update the 'TRUN' column based on dict_consensus (Alu_Seq and L1_Seq lengths)
    df_insertions = update_trun(df_insertions, dict_consensus)

    # Remove possible rows where the 'Length' column is negative
    df_insertions = df_insertions[df_insertions['Length'] >= 0]

    return df_insertions

def add_source_gene_info(df_insertions, source_L1_path, source_SVA_path):
    # Load the source element tables
    table_source_L1 = pd.read_csv(source_L1_path, sep='\t')
    table_source_SVA = pd.read_csv(source_SVA_path, sep='\t')

    # Add necessary columns with default values
    df_insertions[['SRC_identifier', 'SRC_ref', 'SRC_beg', 'SRC_end', 'SRC_cont_PCAWG', 'SRC_strand', 'SRC_in_ref_genome']] = 0

    # Iterate over rows to assign values from table_source_L1 and table_source_SVA
    for index, row in df_insertions.iterrows():
        # If 'L1' AND 'TD' in name, assign values from L1 table
        if 'L1' in row['name'] and 'TD' in row['name']:
            probabilities = table_source_L1['SRC_cont_PCAWG'].values
            selected_row = table_source_L1.sample(weights=probabilities).iloc[0]
            df_insertions.loc[index, ['SRC_identifier', 'SRC_ref', 'SRC_beg', 'SRC_end', 'SRC_cont_PCAWG', 'SRC_strand', 'SRC_in_ref_genome']] = selected_row.values
        # If 'orphan' AND 'TD' in name, assign values assuming orphan gets from L1
        elif row['name'] == 'orphan':
            probabilities = table_source_L1['SRC_cont_PCAWG'].values  
            selected_row = table_source_L1.sample(weights=probabilities).iloc[0]
            df_insertions.loc[index, ['SRC_identifier', 'SRC_ref', 'SRC_beg', 'SRC_end', 'SRC_cont_PCAWG', 'SRC_strand', 'SRC_in_ref_genome']] = selected_row.values
        # If 'SVA' AND 'TD' in name, assign values from SVA table
        elif 'SVA' in row['name'] and 'TD' in row['name']:
            probabilities = table_source_SVA['SRC_contribution'].values
            selected_row = table_source_SVA.sample(weights=probabilities).iloc[0]
            df_insertions.loc[index, ['SRC_identifier', 'SRC_ref', 'SRC_beg', 'SRC_end', 'SRC_contribution', 'SRC_strand', 'SRC_in_ref_genome']] = selected_row.values

    # Fill NaN values with 0
    df_insertions.fillna(0, inplace=True)

    # Add new columns for transduction start and end
    df_insertions[['TD_beg', 'TD_end']] = 0

    # Modify the values of 'TD_beg' and 'TD_end' based on conditions
    for index, row in df_insertions.iterrows():
        if 'L1' in row['name'] and 'TD' in row['name']:
            if row['SRC_strand'] == 'plus':
                df_insertions.at[index, 'TD_beg'] = row['SRC_end']
                df_insertions.at[index, 'TD_end'] = row['SRC_end'] + row['TD_5'] if row['TD_5'] != 0 else row['SRC_end'] + row['TD_3']
            elif row['SRC_strand'] == 'minus':
                df_insertions.at[index, 'TD_end'] = row['SRC_beg']
                df_insertions.at[index, 'TD_beg'] = row['SRC_beg'] - row['TD_5'] if row['TD_5'] != 0 else row['SRC_beg'] - row['TD_3']
        elif 'SVA' in row['name'] and 'TD' in row['name']:
            if row['SRC_strand'] == 'plus':
                df_insertions.at[index, 'TD_beg'] = row['SRC_end']
                df_insertions.at[index, 'TD_end'] = row['SRC_end'] + row['TD_5'] if row['TD_5'] != 0 else row['SRC_end'] + row['TD_3']
            elif row['SRC_strand'] == 'minus':
                df_insertions.at[index, 'TD_end'] = row['SRC_beg']
                df_insertions.at[index, 'TD_beg'] = row['SRC_beg'] - row['TD_5'] if row['TD_5'] != 0 else row['SRC_beg'] - row['TD_3']
        elif 'orphan' in row['name']:
            if row['SRC_strand'] == 'plus':
                df_insertions.at[index, 'TD_beg'] = row['SRC_end']
                df_insertions.at[index, 'TD_end'] = row['SRC_end'] + row['TD_orphan_Length']
            elif row['SRC_strand'] == 'minus':
                df_insertions.at[index, 'TD_end'] = row['SRC_beg']
                df_insertions.at[index, 'TD_beg'] = row['SRC_beg'] - row['TD_orphan_Length']

    # Update SRC_strand values
    df_insertions['SRC_strand'] = df_insertions['SRC_strand'].replace({'minus': '-', 'plus': '+'})

    # Convert 0s to 'NA', except for 'SRC_in_ref_genome' column
    keep_columns = ['SRC_in_ref_genome']
    for column in df_insertions.columns:
        if column not in keep_columns:
            df_insertions[column] = df_insertions[column].apply(lambda x: 'NA' if x == 0 else x)

    # Directly update 'SRC_in_ref_genome' based on 'SRC_identifier'
    df_insertions['SRC_in_ref_genome'] = df_insertions.apply(lambda row: 'NA' if row['SRC_identifier'] == 'NA' else row['SRC_in_ref_genome'], axis=1)

    return df_insertions

# Open and get SVAs VNTR motifs
def read_file_and_store_lines(file_path):
    with open(file_path, 'r') as file:
        lines = [line.strip() for line in file]
    return lines
def VNTR_insertions(row, motifs_file):
    '''
    Function to generate VNTR sequences, selecting motifs based on their proportions
    from a .tsv file containing 'Motif' and 'Proportion' columns.
    '''
    # Step 1: Read the motifs and their proportions from the .tsv file
    motifs_df_original = pd.read_csv(motifs_file, sep='\t')  # Reads .tsv file
    motifs_df = motifs_df_original.copy()  # Create a copy of the DataFrame to work with
    
    # Create a list of motifs and their associated proportions
    motifs = motifs_df['Motif'].tolist()
    proportions = motifs_df['Proportion'].tolist()
    
    # Step 2: Extract the necessary values from the DataFrame row
    total_length = int(row['Length'])
    number_motifs = int(row['VNTR_Num_Motifs'])
    
    selected_motifs = []
    total_selected_length = 0
    
    # Step 3: Select motifs until we have the desired number of motifs and the total length is not exceeded
    while len(selected_motifs) < number_motifs:
        # Randomly select a motif based on the proportions
        motif_index = random.choices(range(len(motifs)), proportions)[0]  # Get the index of the selected motif
        motif = motifs[motif_index]  # Use the index to get the motif
        motif_length = len(motif)
        
        # Check if adding this motif would exceed the total length
        if total_selected_length + motif_length <= total_length:
            selected_motifs.append(motif)
            total_selected_length += motif_length
            
            # Remove the selected motif and its proportion from the DataFrame
            motifs_df = motifs_df.drop(index=motif_index).reset_index(drop=True)  # Remove the selected motif row
            motifs = motifs_df['Motif'].tolist()  # Update the motifs list
            proportions = motifs_df['Proportion'].tolist()  # Update the proportions list
    
    # Step 4: Handle the case of a single motif
    if number_motifs == 1:
        single_motif = selected_motifs[0]
        single_motif_length = len(single_motif)
        repetitions = total_length // single_motif_length
        remainder = total_length % single_motif_length
        
        sequence = single_motif * repetitions
        if remainder > 0:
            sequence += single_motif[:remainder]  # Add partial motif
    
    else:
        # Step 5: If more than one motif, divide total length into non-equal parts
        lengths = random.sample(range(1, total_length // number_motifs + 1), number_motifs - 1)
        lengths.append(total_length - sum(lengths))  # Make sure total length is preserved
        
        lengths = [length for length in lengths if length > 0]  # Filter out any zero-length
        
        sequence = ''
        
        for i in range(number_motifs):
            motif = selected_motifs[i]
            motif_length = len(motif)
            part_length = lengths[i]
            repetitions = part_length // motif_length
            remainder = part_length % motif_length
            
            sequence += motif * repetitions
            if remainder > 0:
                sequence += motif[:remainder]  # Add partial motif
    
    return sequence

def DUP_insertions(row,reference_fasta):
    ''' 
    Function to generate the duplicated sequences
    '''
    start = int(row['beg'])
    length = int(row['Length'])
    end = start + length
    with pysam.FastaFile(reference_fasta) as fasta_file:
        insertion = fasta_file.fetch(row['#ref'], start, end)
    return insertion

def NUMT_insertions(row,dict_consensus):
    ''' 
    Function to generate the mitochondrial insertion sequences
    '''
    # Retrieve the sequence from the dictionary
    seq = dict_consensus['NUMT_Seq']
    
    # Transform Total_Length to integer
    length = int(row['Length'])
    
    # Get a random starting position from the sequence
    start_pos = random.randint(0, len(seq) - 1) #if seq else 0
    
    # Take 'length' number of positions from 'seq', wrapping around if necessary
    result = ''
    for offset in range(length):
        result += seq[(start_pos + offset) % len(seq)]
    
    return result

def inverse_sequence(seq):
    ''' 
    Function to reverse the order of the given sequence
    '''
    # Reverse the sequence using slicing [::-1]
    reversed_seq = seq[::-1]
    
    return reversed_seq

def orphan_insertions(row, reference_fasta):
    ''' 
    Function to generate orphan sequences
    '''
    # Poly A 1
    if pd.isna(row['PolyA_Length_1']) or row['PolyA_Length_1'] == 'NA':
        polyA_length1 = 0 
    else:
        polyA_length1 = int(row['PolyA_Length_1'])
    polyA1 = 'A' * polyA_length1

    # Transduction
    # Fetch a sequence using pysam
    with pysam.FastaFile(reference_fasta) as fasta_file:
        transduction = fasta_file.fetch(row['SRC_ref'], row['TD_beg'], row['TD_end'])

    # TSD
    beg = int(row['beg'])
    tsd_len = int(row['TSD_Length'])
    start_tsd = beg - tsd_len
    # Fetch a sequence using pysam
    with pysam.FastaFile(reference_fasta) as fasta_file:
        TSD = fasta_file.fetch(row['#ref'], start_tsd, beg)

    # Final sequence
    seq = transduction + polyA1 + TSD
    
    return seq

def Alu__FOR_POLYA_insertions(row, dict_consensus, reference_fasta):
    # Extract the value from the 'FOR'
    for_value = row['FOR']
    # Get the Alu_Seq from the dictionary
    Alu_consensus = dict_consensus.get('Alu_Seq', '')
    
    # Slice the Alu_Seq from the end based on the 'FOR' value
    # Ensure that for_value is not greater than the length of Alu_Seq
    if isinstance(for_value, int) and for_value <= len(Alu_consensus):
        Alu_seq = Alu_consensus[-for_value:]
    else:
        # If for_value is larger than the length of the sequence, return the whole sequence
        Alu_seq = Alu_consensus
    
    # Poly A 1
    if pd.isna(row['PolyA_Length_1']) or row['PolyA_Length_1'] == 'NA':
        polyA_length1 = 0 
    else:
        polyA_length1 = int(row['PolyA_Length_1'])
    polyA1 = 'A' * polyA_length1

    # TSD
    beg = int(row['beg'])
    tsd_len = int(row['TSD_Length'])
    start_tsd = beg - tsd_len
    # Fetch a sequence using pysam
    with pysam.FastaFile(reference_fasta) as fasta_file:
        TSD = fasta_file.fetch(row['#ref'], start_tsd, beg)

    # Final sequence
    seq = Alu_seq + polyA1 + TSD

    return seq

def L1__FOR_POLYA_insertions(row, dict_consensus, reference_fasta):
    # Extract the value from the 'FOR'
    for_value = row['FOR']
    # Get the L1_Seq from the dictionary
    L1_consensus = dict_consensus.get('L1_Seq', '')
    
    # Slice the L1_Seq from the end based on the 'FOR' value
    # Ensure that for_value is not greater than the length of L1_Seq
    if isinstance(for_value, int) and for_value <= len(L1_consensus):
        L1_seq = L1_consensus[-for_value:]
    else:
        # If for_value is larger than the length of the sequence, return the whole sequence
        L1_seq = L1_consensus
    
    # Poly A 1
    if pd.isna(row['PolyA_Length_1']) or row['PolyA_Length_1'] == 'NA':
        polyA_length1 = 0 
    else:
        polyA_length1 = int(row['PolyA_Length_1'])
    polyA1 = 'A' * polyA_length1

    # TSD
    beg = int(row['beg'])
    tsd_len = int(row['TSD_Length'])
    start_tsd = beg - tsd_len
    # Fetch a sequence using pysam
    with pysam.FastaFile(reference_fasta) as fasta_file:
        TSD = fasta_file.fetch(row['#ref'], start_tsd, beg)

    # Final sequence
    seq = L1_seq + polyA1 + TSD

    return seq

def L1__TD_FOR_POLYA_insertions(row, dict_consensus, reference_fasta):
    # Extract the value from the 'FOR'
    for_value = row['FOR']
    # Get the L1_Seq from the dictionary
    L1_consensus = dict_consensus.get('L1_Seq', '')
    
    # Slice the L1_Seq from the end based on the 'FOR' value
    # Ensure that for_value is not greater than the length of L1_Seq
    if isinstance(for_value, int) and for_value <= len(L1_consensus):
        L1_seq = L1_consensus[-for_value:]
    else:
        # If for_value is larger than the length of the sequence, return the whole sequence
        L1_seq = L1_consensus
    
    # Poly A 1
    if pd.isna(row['PolyA_Length_1']) or row['PolyA_Length_1'] == 'NA':
        polyA_length1 = 0 
    else:
        polyA_length1 = int(row['PolyA_Length_1'])
    polyA1 = 'A' * polyA_length1

    # TSD
    beg = int(row['beg'])
    tsd_len = int(row['TSD_Length'])
    start_tsd = beg - tsd_len
    # Fetch a sequence using pysam
    with pysam.FastaFile(reference_fasta) as fasta_file:
        TSD = fasta_file.fetch(row['#ref'], start_tsd, beg)

    # Transduction
    # Fetch a sequence using pysam
    with pysam.FastaFile(reference_fasta) as fasta_file:
        transduction = fasta_file.fetch(row['SRC_ref'], row['TD_beg'], row['TD_end'])

    # Final sequence
    seq = transduction + L1_seq + polyA1 + TSD

    return seq

def L1__FOR_POLYA_TD_POLYA_insertions(row, dict_consensus, reference_fasta):
    # Extract the value from the 'FOR'
    for_value = row['FOR']
    # Get the L1_Seq from the dictionary
    L1_consensus = dict_consensus.get('L1_Seq', '')
    
    # Slice the L1_Seq from the end based on the 'FOR' value
    # Ensure that for_value is not greater than the length of L1_Seq
    if isinstance(for_value, int) and for_value <= len(L1_consensus):
        L1_seq = L1_consensus[-for_value:]
    else:
        # If for_value is larger than the length of the sequence, return the whole sequence
        L1_seq = L1_consensus
    
    # Poly A 1
    if pd.isna(row['PolyA_Length_1']) or row['PolyA_Length_1'] == 'NA':
        polyA_length1 = 0 
    else:
        polyA_length1 = int(row['PolyA_Length_1'])
    polyA1 = 'A' * polyA_length1

    # Poly A 2
    if pd.isna(row['PolyA_Length_2']) or row['PolyA_Length_2'] == 'NA':
        polyA_length2 = 0 
    else:
        polyA_length2 = int(row['PolyA_Length_2'])
    polyA2 = 'A' * polyA_length2

    # TSD
    beg = int(row['beg'])
    tsd_len = int(row['TSD_Length'])
    start_tsd = beg - tsd_len
    # Fetch a sequence using pysam
    with pysam.FastaFile(reference_fasta) as fasta_file:
        TSD = fasta_file.fetch(row['#ref'], start_tsd, beg)

    # Transduction
    # Fetch a sequence using pysam
    with pysam.FastaFile(reference_fasta) as fasta_file:
        transduction = fasta_file.fetch(row['SRC_ref'], row['TD_beg'], row['TD_end'])
        
    # Final sequence
    seq = L1_seq + polyA1 + transduction + polyA2 + TSD

    return seq

def L1__TRUN_REV_BLUNT_FOR_POLYA_insertions(row, dict_consensus, reference_fasta):
    # Extract the value from the 'FOR' and 'REV' columns
    for_value = row['FOR']
    rev_value = row['REV'] 
    
    # Get the L1_Seq from the dictionary
    L1_consensus = dict_consensus.get('L1_Seq', '')
    
    # Slice the L1_Seq from the end based on the 'FOR' value
    if isinstance(for_value, int) and for_value <= len(L1_consensus):
        L1_seq = L1_consensus[-for_value:]
    else:
        L1_seq = L1_consensus
    
    # Get the REV_seq based on the 'REV' value
    if isinstance(rev_value, int) and rev_value <= len(L1_consensus):
        # Slice the sequence from the end based on the 'REV' value starting after L1_seq
        rev_seq_start = len(L1_consensus) - for_value  # The point where L1_seq ends
        rev_seq = L1_consensus[rev_seq_start - rev_value: rev_seq_start]
        rev_seq = inverse_sequence(rev_seq)  # Apply inverse function to REV sequence
    else:
        rev_seq = ''  # In case REV value is invalid or too large
    
    # Poly A 1
    if pd.isna(row['PolyA_Length_1']) or row['PolyA_Length_1'] == 'NA':
        polyA_length1 = 0 
    else:
        polyA_length1 = int(row['PolyA_Length_1'])
    polyA1 = 'A' * polyA_length1

    # TSD
    beg = int(row['beg'])
    tsd_len = int(row['TSD_Length'])
    start_tsd = beg - tsd_len
    # Fetch a sequence using pysam
    with pysam.FastaFile(reference_fasta) as fasta_file:
        TSD = fasta_file.fetch(row['#ref'], start_tsd, beg)

    # Final sequence
    seq = rev_seq + L1_seq + polyA1 + TSD
    
    return seq

def L1__TRUN_REV_BLUNT_FOR_POLYA_TD_POLYA_insertions(row, dict_consensus, reference_fasta):
    # Extract the value from the 'FOR' and 'REV' columns
    for_value = row['FOR']
    rev_value = row['REV'] 
    
    # Get the L1_Seq from the dictionary
    L1_consensus = dict_consensus.get('L1_Seq', '')
    
    # Slice the L1_Seq from the end based on the 'FOR' value
    if isinstance(for_value, int) and for_value <= len(L1_consensus):
        L1_seq = L1_consensus[-for_value:]
    else:
        L1_seq = L1_consensus
    
    # Get the REV_seq based on the 'REV' value
    if isinstance(rev_value, int) and rev_value <= len(L1_consensus):
        # Slice the sequence from the end based on the 'REV' value starting after L1_seq
        rev_seq_start = len(L1_consensus) - for_value  # The point where L1_seq ends
        rev_seq = L1_consensus[rev_seq_start - rev_value: rev_seq_start]
        rev_seq = inverse_sequence(rev_seq)  # Apply inverse function to REV sequence
    else:
        rev_seq = ''  # In case REV value is invalid or too large
    
    # Poly A 1
    if pd.isna(row['PolyA_Length_1']) or row['PolyA_Length_1'] == 'NA':
        polyA_length1 = 0 
    else:
        polyA_length1 = int(row['PolyA_Length_1'])
    polyA1 = 'A' * polyA_length1

    # Poly A 2
    if pd.isna(row['PolyA_Length_2']) or row['PolyA_Length_2'] == 'NA':
        polyA_length2 = 0 
    else:
        polyA_length2 = int(row['PolyA_Length_2'])
    polyA2 = 'A' * polyA_length2

    # Transduction
    # Fetch a sequence using pysam
    with pysam.FastaFile(reference_fasta) as fasta_file:
        transduction = fasta_file.fetch(row['SRC_ref'], row['TD_beg'], row['TD_end'])

    # TSD
    beg = int(row['beg'])
    tsd_len = int(row['TSD_Length'])
    start_tsd = beg - tsd_len
    # Fetch a sequence using pysam
    with pysam.FastaFile(reference_fasta) as fasta_file:
        TSD = fasta_file.fetch(row['#ref'], start_tsd, beg)

    # Final sequence
    seq = rev_seq + L1_seq + polyA1 + transduction + polyA2 + TSD
    
    return seq

def L1__TRUN_REV_DUP_FOR_POLYA_insertions(row, dict_consensus, reference_fasta):
    # Extract the value from the 'FOR', 'REV', and 'DUP' columns
    for_value = row['FOR']
    rev_value = row['REV']
    dup_value = row['DUP']
    
    # Get the L1_Seq from the dictionary
    L1_consensus = dict_consensus.get('L1_Seq', '')
    
    # Slice the L1_Seq from the end based on the 'FOR' value
    if isinstance(for_value, int) and for_value <= len(L1_consensus):
        L1_seq = L1_consensus[-for_value:]
    else:
        L1_seq = L1_consensus
        
    # Get the REV_seq based on the 'REV' value
    if isinstance(rev_value, int) and rev_value <= len(L1_consensus):
        # Slice the sequence from the end based on the 'REV' value starting after L1_seq
        rev_seq_start = len(L1_consensus) - for_value  # The point where L1_seq ends
        rev_seq = L1_consensus[rev_seq_start - rev_value: rev_seq_start]
        rev_seq = inverse_sequence(rev_seq)  # Apply inverse function to REV sequence
    else:
        rev_seq = ''  # In case REV value is invalid or too large
    
    # Create DUP_seq by duplicating the 'DUP' bases of L1_seq 
    if isinstance(dup_value, int) and dup_value > 0:
        dup_seq = L1_seq[:dup_value]
    else:
        dup_seq = ''  # If DUP value is invalid or 0
    
    # Poly A 1
    if pd.isna(row['PolyA_Length_1']) or row['PolyA_Length_1'] == 'NA':
        polyA_length1 = 0
    else:
        polyA_length1 = int(row['PolyA_Length_1'])
    polyA1 = 'A' * polyA_length1

    # TSD
    beg = int(row['beg'])
    tsd_len = int(row['TSD_Length'])
    start_tsd = beg - tsd_len
    # Fetch a sequence using pysam
    with pysam.FastaFile(reference_fasta) as fasta_file:
        TSD = fasta_file.fetch(row['#ref'], start_tsd, beg)

    # Final sequence
    seq = rev_seq + dup_seq + L1_seq + polyA1 + TSD
    
    return seq

def L1__TRUN_REV_DUP_FOR_POLYA_TD_POLYA_insertions(row, dict_consensus, reference_fasta):
    # Extract the value from the 'FOR', 'REV', and 'DUP' columns
    for_value = row['FOR']
    rev_value = row['REV']
    dup_value = row['DUP']
    
    # Get the L1_Seq from the dictionary
    L1_consensus = dict_consensus.get('L1_Seq', '')
    
    # Slice the L1_Seq from the end based on the 'FOR' value
    if isinstance(for_value, int) and for_value <= len(L1_consensus):
        L1_seq = L1_consensus[-for_value:]
    else:
        L1_seq = L1_consensus
        
    # Get the REV_seq based on the 'REV' value
    if isinstance(rev_value, int) and rev_value <= len(L1_consensus):
        # Slice the sequence from the end based on the 'REV' value starting after L1_seq
        rev_seq_start = len(L1_consensus) - for_value  # The point where L1_seq ends
        rev_seq = L1_consensus[rev_seq_start - rev_value: rev_seq_start]
        rev_seq = inverse_sequence(rev_seq)  # Apply inverse function to REV sequence
    else:
        rev_seq = ''  # In case REV value is invalid or too large
    
    # Create DUP_seq by duplicating the 'DUP' bases of L1_seq 
    if isinstance(dup_value, int) and dup_value > 0:
        dup_seq = L1_seq[:dup_value]
    else:
        dup_seq = ''  # If DUP value is invalid or 0
    
    # Poly A 1
    if pd.isna(row['PolyA_Length_1']) or row['PolyA_Length_1'] == 'NA':
        polyA_length1 = 0
    else:
        polyA_length1 = int(row['PolyA_Length_1'])
    polyA1 = 'A' * polyA_length1

    # Poly A 2
    if pd.isna(row['PolyA_Length_2']) or row['PolyA_Length_2'] == 'NA':
        polyA_length2 = 0 
    else:
        polyA_length2 = int(row['PolyA_Length_2'])
    polyA2 = 'A' * polyA_length2

    # Transduction
    # Fetch a sequence using pysam
    with pysam.FastaFile(reference_fasta) as fasta_file:
        transduction = fasta_file.fetch(row['SRC_ref'], row['TD_beg'], row['TD_end'])

    # TSD
    beg = int(row['beg'])
    tsd_len = int(row['TSD_Length'])
    start_tsd = beg - tsd_len
    # Fetch a sequence using pysam
    with pysam.FastaFile(reference_fasta) as fasta_file:
        TSD = fasta_file.fetch(row['#ref'], start_tsd, beg)

    # Final sequence
    seq = rev_seq + dup_seq + L1_seq + polyA1 + transduction + polyA2 + TSD
    
    return seq

def L1__TRUN_REV_DEL_FOR_POLYA_insertions(row, dict_consensus, reference_fasta):
    # Extract the value from the 'FOR', 'REV', and 'DUP' columns
    for_value = row['FOR']
    rev_value = row['REV']
    del_value = row['DEL']
    
    # Get the L1_Seq from the dictionary
    L1_consensus = dict_consensus.get('L1_Seq', '')
    
    # Slice the L1_Seq from the end based on 'FOR' - 'DEL'
    if isinstance(for_value, int) and for_value > del_value and (for_value - del_value) <= len(L1_consensus):
        L1_seq = L1_consensus[-(for_value - del_value):]
    else:
        L1_seq = L1_consensus
        
    # Get the REV_seq based on the 'REV' value
    if isinstance(rev_value, int) and rev_value <= len(L1_consensus):
        # Slice the sequence from the end based on the 'REV' value starting after L1_seq
        rev_seq_start = len(L1_consensus) - for_value  # The point where L1_seq ends
        rev_seq = L1_consensus[rev_seq_start - rev_value: rev_seq_start]
        rev_seq = inverse_sequence(rev_seq)  # Apply inverse function to REV sequence
    else:
        rev_seq = ''  # In case REV value is invalid or too large
    
    # Poly A 1
    if pd.isna(row['PolyA_Length_1']) or row['PolyA_Length_1'] == 'NA':
        polyA_length1 = 0
    else:
        polyA_length1 = int(row['PolyA_Length_1'])
    polyA1 = 'A' * polyA_length1

    # TSD
    beg = int(row['beg'])
    tsd_len = int(row['TSD_Length'])
    start_tsd = beg - tsd_len
    # Fetch a sequence using pysam
    with pysam.FastaFile(reference_fasta) as fasta_file:
        TSD = fasta_file.fetch(row['#ref'], start_tsd, beg)

    # Final sequence
    seq = rev_seq + L1_seq + polyA1 + TSD
    
    return seq

def L1__TRUN_REV_DEL_FOR_POLYA_TD_POLYA_insertions(row, dict_consensus, reference_fasta):
    # Extract the value from the 'FOR', 'REV', and 'DUP' columns
    for_value = row['FOR']
    rev_value = row['REV']
    del_value = row['DEL']
    
    # Get the L1_Seq from the dictionary
    L1_consensus = dict_consensus.get('L1_Seq', '')
    
    # Slice the L1_Seq from the end based on 'FOR' - 'DEL'
    if isinstance(for_value, int) and for_value > del_value and (for_value - del_value) <= len(L1_consensus):
        L1_seq = L1_consensus[-(for_value - del_value):]
    else:
        L1_seq = L1_consensus
        
    # Get the REV_seq based on the 'REV' value
    if isinstance(rev_value, int) and rev_value <= len(L1_consensus):
        # Slice the sequence from the end based on the 'REV' value starting after L1_seq
        rev_seq_start = len(L1_consensus) - for_value  # The point where L1_seq ends
        rev_seq = L1_consensus[rev_seq_start - rev_value: rev_seq_start]
        rev_seq = inverse_sequence(rev_seq)  # Apply inverse function to REV sequence
    else:
        rev_seq = ''  # In case REV value is invalid or too large
    
    # Poly A 1
    if pd.isna(row['PolyA_Length_1']) or row['PolyA_Length_1'] == 'NA':
        polyA_length1 = 0
    else:
        polyA_length1 = int(row['PolyA_Length_1'])
    polyA1 = 'A' * polyA_length1

    # Poly A 2
    if pd.isna(row['PolyA_Length_2']) or row['PolyA_Length_2'] == 'NA':
        polyA_length2 = 0 
    else:
        polyA_length2 = int(row['PolyA_Length_2'])
    polyA2 = 'A' * polyA_length2

    # Transduction
    # Fetch a sequence using pysam
    with pysam.FastaFile(reference_fasta) as fasta_file:
        transduction = fasta_file.fetch(row['SRC_ref'], row['TD_beg'], row['TD_end'])

    # TSD
    beg = int(row['beg'])
    tsd_len = int(row['TSD_Length'])
    start_tsd = beg - tsd_len
    # Fetch a sequence using pysam
    with pysam.FastaFile(reference_fasta) as fasta_file:
        TSD = fasta_file.fetch(row['#ref'], start_tsd, beg)

    # Final sequence
    seq = rev_seq + L1_seq + polyA1 + transduction + polyA2 + TSD
    
    return seq

def SVA__SINE_R_POLYA_insertions(row, dict_consensus, reference_fasta):
    # Get the L1_Seq from the dictionary
    SINE_R_seq = dict_consensus.get('SVA_SINE-R_Seq', '')
    
    # Poly A 1
    if pd.isna(row['PolyA_Length_1']) or row['PolyA_Length_1'] == 'NA':
        polyA_length1 = 0
    else:
        polyA_length1 = int(row['PolyA_Length_1'])
    polyA1 = 'A' * polyA_length1

    # TSD
    beg = int(row['beg'])
    tsd_len = int(row['TSD_Length'])
    start_tsd = beg - tsd_len
    # Fetch a sequence using pysam
    with pysam.FastaFile(reference_fasta) as fasta_file:
        TSD = fasta_file.fetch(row['#ref'], start_tsd, beg)

    # Final sequence
    seq = SINE_R_seq + polyA1 + TSD
    
    return seq

def SVA__VNTR_SINE_R_POLYA_insertions(row, dict_consensus, reference_fasta, SVA_VNTR_list):
    # Get the L1_Seq from the dictionary
    SINE_R_seq = dict_consensus.get('SVA_SINE-R_Seq', '')
    
    # Poly A 1
    if pd.isna(row['PolyA_Length_1']) or row['PolyA_Length_1'] == 'NA':
        polyA_length1 = 0
    else:
        polyA_length1 = int(row['PolyA_Length_1'])
    polyA1 = 'A' * polyA_length1

    # TSD
    beg = int(row['beg'])
    tsd_len = int(row['TSD_Length'])
    start_tsd = beg - tsd_len
    # Fetch a sequence using pysam
    with pysam.FastaFile(reference_fasta) as fasta_file:
        TSD = fasta_file.fetch(row['#ref'], start_tsd, beg)

    # VNTR
    random_row = random.choice(SVA_VNTR_list).strip()  # Select random VNTR motif
    sva_vntr_length = int(row['SVA_VNTR_Length'])  # Desired length of VNTR
    
    # Generate VNTR sequence taking desired bases
    # In case VNTR is smaller than desired bases, repeat and wrap around the selected sequence
    vntr_sequence = ''
    while len(vntr_sequence) < sva_vntr_length:
        vntr_sequence += random_row
    vntr_sequence = vntr_sequence[:sva_vntr_length]

    # Final sequence
    seq = vntr_sequence + SINE_R_seq + polyA1 + TSD
    
    return seq

def SVA__Alu_like_VNTR_SINE_R_POLYA_insertions(row, dict_consensus, reference_fasta, SVA_VNTR_list):
    # Get the L1_Seq from the dictionary
    SINE_R_seq = dict_consensus.get('SVA_SINE-R_Seq', '')
    Alu_like_seq = dict_consensus.get('SVA_Alu-like_Seq', '')

    # Poly A 1
    if pd.isna(row['PolyA_Length_1']) or row['PolyA_Length_1'] == 'NA':
        polyA_length1 = 0
    else:
        polyA_length1 = int(row['PolyA_Length_1'])
    polyA1 = 'A' * polyA_length1

    # TSD
    beg = int(row['beg'])
    tsd_len = int(row['TSD_Length'])
    start_tsd = beg - tsd_len
    # Fetch a sequence using pysam
    with pysam.FastaFile(reference_fasta) as fasta_file:
        TSD = fasta_file.fetch(row['#ref'], start_tsd, beg)

    # VNTR
    random_row = random.choice(SVA_VNTR_list).strip()  # Select random VNTR motif
    sva_vntr_length = int(row['SVA_VNTR_Length'])  # Desired length of VNTR
    
    # Generate VNTR sequence taking desired bases
    # In case VNTR is smaller than desired bases, repeat and wrap around the selected sequence
    vntr_sequence = ''
    while len(vntr_sequence) < sva_vntr_length:
        vntr_sequence += random_row
    vntr_sequence = vntr_sequence[:sva_vntr_length]

    # Final sequence
    seq = Alu_like_seq + vntr_sequence + SINE_R_seq + polyA1 + TSD
    
    return seq

def SVA__MAST2_VNTR_SINE_R_POLYA_insertions(row, dict_consensus, reference_fasta, SVA_VNTR_list):
    # Get the L1_Seq from the dictionary
    SINE_R_seq = dict_consensus.get('SVA_SINE-R_Seq', '')
    MAST2_seq = dict_consensus.get('SVA_MAST2_Seq', '')

    # Poly A 1
    if pd.isna(row['PolyA_Length_1']) or row['PolyA_Length_1'] == 'NA':
        polyA_length1 = 0
    else:
        polyA_length1 = int(row['PolyA_Length_1'])
    polyA1 = 'A' * polyA_length1

    # TSD
    beg = int(row['beg'])
    tsd_len = int(row['TSD_Length'])
    start_tsd = beg - tsd_len
    # Fetch a sequence using pysam
    with pysam.FastaFile(reference_fasta) as fasta_file:
        TSD = fasta_file.fetch(row['#ref'], start_tsd, beg)

    # VNTR
    random_row = random.choice(SVA_VNTR_list).strip()  # Select random VNTR motif
    sva_vntr_length = int(row['SVA_VNTR_Length'])  # Desired length of VNTR
    
    # Generate VNTR sequence taking desired bases
    # In case VNTR is smaller than desired bases, repeat and wrap around the selected sequence
    vntr_sequence = ''
    while len(vntr_sequence) < sva_vntr_length:
        vntr_sequence += random_row
    vntr_sequence = vntr_sequence[:sva_vntr_length]

    # Final sequence
    seq = MAST2_seq + vntr_sequence + SINE_R_seq + polyA1 + TSD
    
    return seq

def SVA__TD_MAST2_VNTR_SINE_R_POLYA_insertions(row, dict_consensus, reference_fasta, SVA_VNTR_list):
    # Get the L1_Seq from the dictionary
    SINE_R_seq = dict_consensus.get('SVA_SINE-R_Seq', '')
    MAST2_seq = dict_consensus.get('SVA_MAST2_Seq', '')

    # Poly A 1
    if pd.isna(row['PolyA_Length_1']) or row['PolyA_Length_1'] == 'NA':
        polyA_length1 = 0
    else:
        polyA_length1 = int(row['PolyA_Length_1'])
    polyA1 = 'A' * polyA_length1

    # TSD
    beg = int(row['beg'])
    tsd_len = int(row['TSD_Length'])
    start_tsd = beg - tsd_len
    # Fetch a sequence using pysam
    with pysam.FastaFile(reference_fasta) as fasta_file:
        TSD = fasta_file.fetch(row['#ref'], start_tsd, beg)

    # VNTR
    random_row = random.choice(SVA_VNTR_list).strip()  # Select random VNTR motif
    sva_vntr_length = int(row['SVA_VNTR_Length'])  # Desired length of VNTR
    
    # Generate VNTR sequence taking desired bases
    # In case VNTR is smaller than desired bases, repeat and wrap around the selected sequence
    vntr_sequence = ''
    while len(vntr_sequence) < sva_vntr_length:
        vntr_sequence += random_row
    vntr_sequence = vntr_sequence[:sva_vntr_length]

    # Transduction
    # Fetch a sequence using pysam
    with pysam.FastaFile(reference_fasta) as fasta_file:
        transduction = fasta_file.fetch(row['SRC_ref'], row['TD_beg'], row['TD_end'])

    # Final sequence
    seq = transduction + MAST2_seq + vntr_sequence + SINE_R_seq + polyA1 + TSD
    
    return seq

def SVA__SINE_R_POLYA_TD_POLYA_insertions(row, dict_consensus, reference_fasta):
    # Get the L1_Seq from the dictionary
    SINE_R_seq = dict_consensus.get('SVA_SINE-R_Seq', '')

    # Poly A 1
    if pd.isna(row['PolyA_Length_1']) or row['PolyA_Length_1'] == 'NA':
        polyA_length1 = 0
    else:
        polyA_length1 = int(row['PolyA_Length_1'])
    polyA1 = 'A' * polyA_length1

    # Poly A 2
    if pd.isna(row['PolyA_Length_2']) or row['PolyA_Length_2'] == 'NA':
        polyA_length2 = 0 
    else:
        polyA_length2 = int(row['PolyA_Length_2'])
    polyA2 = 'A' * polyA_length2

    # TSD
    beg = int(row['beg'])
    tsd_len = int(row['TSD_Length'])
    start_tsd = beg - tsd_len
    # Fetch a sequence using pysam
    with pysam.FastaFile(reference_fasta) as fasta_file:
        TSD = fasta_file.fetch(row['#ref'], start_tsd, beg)


    # Transduction
    # Fetch a sequence using pysam
    with pysam.FastaFile(reference_fasta) as fasta_file:
        transduction = fasta_file.fetch(row['SRC_ref'], row['TD_beg'], row['TD_end'])

    # Final sequence
    seq = SINE_R_seq + polyA1 + transduction + polyA2 + TSD

    return seq

def SVA__VNTR_SINE_R_POLYA_TD_POLYA_insertions(row, dict_consensus, reference_fasta, SVA_VNTR_list):
    # Get the L1_Seq from the dictionary
    SINE_R_seq = dict_consensus.get('SVA_SINE-R_Seq', '')

    # Poly A 1
    if pd.isna(row['PolyA_Length_1']) or row['PolyA_Length_1'] == 'NA':
        polyA_length1 = 0
    else:
        polyA_length1 = int(row['PolyA_Length_1'])
    polyA1 = 'A' * polyA_length1

    # Poly A 2
    if pd.isna(row['PolyA_Length_2']) or row['PolyA_Length_2'] == 'NA':
        polyA_length2 = 0 
    else:
        polyA_length2 = int(row['PolyA_Length_2'])
    polyA2 = 'A' * polyA_length2

    # TSD
    beg = int(row['beg'])
    tsd_len = int(row['TSD_Length'])
    start_tsd = beg - tsd_len
    # Fetch a sequence using pysam
    with pysam.FastaFile(reference_fasta) as fasta_file:
        TSD = fasta_file.fetch(row['#ref'], start_tsd, beg)

    # VNTR
    random_row = random.choice(SVA_VNTR_list).strip()  # Select random VNTR motif
    sva_vntr_length = int(row['SVA_VNTR_Length'])  # Desired length of VNTR
    
    # Generate VNTR sequence taking desired bases
    # In case VNTR is smaller than desired bases, repeat and wrap around the selected sequence
    vntr_sequence = ''
    while len(vntr_sequence) < sva_vntr_length:
        vntr_sequence += random_row
    vntr_sequence = vntr_sequence[:sva_vntr_length]

    # Transduction
    # Fetch a sequence using pysam
    with pysam.FastaFile(reference_fasta) as fasta_file:
        transduction = fasta_file.fetch(row['SRC_ref'], row['TD_beg'], row['TD_end'])

    # Final sequence
    seq = vntr_sequence + SINE_R_seq + polyA1 + transduction + polyA2 + TSD

    return seq

def SVA__Alu_like_VNTR_SINE_R_POLYA_TD_POLYA_insertions(row, dict_consensus, reference_fasta, SVA_VNTR_list):
    # Get the L1_Seq from the dictionary
    SINE_R_seq = dict_consensus.get('SVA_SINE-R_Seq', '')
    Alu_like_seq = dict_consensus.get('SVA_Alu-like_Seq', '')

    # Poly A 1
    if pd.isna(row['PolyA_Length_1']) or row['PolyA_Length_1'] == 'NA':
        polyA_length1 = 0
    else:
        polyA_length1 = int(row['PolyA_Length_1'])
    polyA1 = 'A' * polyA_length1

    # Poly A 2
    if pd.isna(row['PolyA_Length_2']) or row['PolyA_Length_2'] == 'NA':
        polyA_length2 = 0 
    else:
        polyA_length2 = int(row['PolyA_Length_2'])
    polyA2 = 'A' * polyA_length2

    # TSD
    beg = int(row['beg'])
    tsd_len = int(row['TSD_Length'])
    start_tsd = beg - tsd_len
    # Fetch a sequence using pysam
    with pysam.FastaFile(reference_fasta) as fasta_file:
        TSD = fasta_file.fetch(row['#ref'], start_tsd, beg)

    # VNTR
    random_row = random.choice(SVA_VNTR_list).strip()  # Select random VNTR motif
    sva_vntr_length = int(row['SVA_VNTR_Length'])  # Desired length of VNTR
    
    # Generate VNTR sequence taking desired bases
    # In case VNTR is smaller than desired bases, repeat and wrap around the selected sequence
    vntr_sequence = ''
    while len(vntr_sequence) < sva_vntr_length:
        vntr_sequence += random_row
    vntr_sequence = vntr_sequence[:sva_vntr_length]

    # Transduction
    # Fetch a sequence using pysam
    with pysam.FastaFile(reference_fasta) as fasta_file:
        transduction = fasta_file.fetch(row['SRC_ref'], row['TD_beg'], row['TD_end'])

    # Final sequence
    seq = Alu_like_seq + vntr_sequence + SINE_R_seq + polyA1 + transduction + polyA2 + TSD

    return seq

def SVA__MAST2_VNTR_SINE_R_POLYA_TD_POLYA_insertions(row, dict_consensus, reference_fasta, SVA_VNTR_list):
    # Get the L1_Seq from the dictionary
    SINE_R_seq = dict_consensus.get('SVA_SINE-R_Seq', '')
    MAST2_seq = dict_consensus.get('SVA_MAST2_Seq', '')

    # Poly A 1
    if pd.isna(row['PolyA_Length_1']) or row['PolyA_Length_1'] == 'NA':
        polyA_length1 = 0
    else:
        polyA_length1 = int(row['PolyA_Length_1'])
    polyA1 = 'A' * polyA_length1

    # Poly A 2
    if pd.isna(row['PolyA_Length_2']) or row['PolyA_Length_2'] == 'NA':
        polyA_length2 = 0 
    else:
        polyA_length2 = int(row['PolyA_Length_2'])
    polyA2 = 'A' * polyA_length2

    # TSD
    beg = int(row['beg'])
    tsd_len = int(row['TSD_Length'])
    start_tsd = beg - tsd_len
    # Fetch a sequence using pysam
    with pysam.FastaFile(reference_fasta) as fasta_file:
        TSD = fasta_file.fetch(row['#ref'], start_tsd, beg)

    # VNTR
    random_row = random.choice(SVA_VNTR_list).strip()  # Select random VNTR motif
    sva_vntr_length = int(row['SVA_VNTR_Length'])  # Desired length of VNTR
    
    # Generate VNTR sequence taking desired bases
    # In case VNTR is smaller than desired bases, repeat and wrap around the selected sequence
    vntr_sequence = ''
    while len(vntr_sequence) < sva_vntr_length:
        vntr_sequence += random_row
    vntr_sequence = vntr_sequence[:sva_vntr_length]

    # Transduction
    # Fetch a sequence using pysam
    with pysam.FastaFile(reference_fasta) as fasta_file:
        transduction = fasta_file.fetch(row['SRC_ref'], row['TD_beg'], row['TD_end'])

    # Final sequence
    seq = MAST2_seq + vntr_sequence + SINE_R_seq + polyA1 + transduction + polyA2 + TSD

    return seq

def SVA__Hexamer_Alu_like_VNTR_SINE_R_POLYA_insertions(row, dict_consensus, reference_fasta, SVA_VNTR_list):
    # Get the L1_Seq from the dictionary
    SINE_R_seq = dict_consensus.get('SVA_SINE-R_Seq', '')
    Alu_like_seq = dict_consensus.get('SVA_Alu-like_Seq', '')

    # Poly A 1
    if pd.isna(row['PolyA_Length_1']) or row['PolyA_Length_1'] == 'NA':
        polyA_length1 = 0
    else:
        polyA_length1 = int(row['PolyA_Length_1'])
    polyA1 = 'A' * polyA_length1

    # TSD
    beg = int(row['beg'])
    tsd_len = int(row['TSD_Length'])
    start_tsd = beg - tsd_len
    # Fetch a sequence using pysam
    with pysam.FastaFile(reference_fasta) as fasta_file:
        TSD = fasta_file.fetch(row['#ref'], start_tsd, beg)

    # VNTR
    random_row = random.choice(SVA_VNTR_list).strip()  # Select random VNTR motif
    sva_vntr_length = int(row['SVA_VNTR_Length'])  # Desired length of VNTR
    
    # Generate VNTR sequence taking desired bases
    # In case VNTR is smaller than desired bases, repeat and wrap around the selected sequence
    vntr_sequence = ''
    while len(vntr_sequence) < sva_vntr_length:
        vntr_sequence += random_row
    vntr_sequence = vntr_sequence[:sva_vntr_length]

    # Hexamer
    hexamer = 'CCCTCT'
    sva_hexamer_length = int(row['SVA_Hexamer'])  # Desired length of hexamer
    # Repeat the hexamer until we exceed or match the required length, then slice to get the exact length
    hexamer_seq = hexamer * (sva_hexamer_length // 6) + hexamer[:sva_hexamer_length % 6]

    # Final sequence
    seq = hexamer_seq + Alu_like_seq + vntr_sequence + SINE_R_seq + polyA1 + TSD

    return seq

def SVA__TD_Hexamer_Alu_like_VNTR_SINE_R_POLYA_insertions(row, dict_consensus, reference_fasta, SVA_VNTR_list):
    # Get the L1_Seq from the dictionary
    SINE_R_seq = dict_consensus.get('SVA_SINE-R_Seq', '')
    Alu_like_seq = dict_consensus.get('SVA_Alu-like_Seq', '')

    # Poly A 1
    if pd.isna(row['PolyA_Length_1']) or row['PolyA_Length_1'] == 'NA':
        polyA_length1 = 0
    else:
        polyA_length1 = int(row['PolyA_Length_1'])
    polyA1 = 'A' * polyA_length1

    # TSD
    beg = int(row['beg'])
    tsd_len = int(row['TSD_Length'])
    start_tsd = beg - tsd_len
    # Fetch a sequence using pysam
    with pysam.FastaFile(reference_fasta) as fasta_file:
        TSD = fasta_file.fetch(row['#ref'], start_tsd, beg)

    # VNTR
    random_row = random.choice(SVA_VNTR_list).strip()  # Select random VNTR motif
    sva_vntr_length = int(row['SVA_VNTR_Length'])  # Desired length of VNTR
    
    # Generate VNTR sequence taking desired bases
    # In case VNTR is smaller than desired bases, repeat and wrap around the selected sequence
    vntr_sequence = ''
    while len(vntr_sequence) < sva_vntr_length:
        vntr_sequence += random_row
    vntr_sequence = vntr_sequence[:sva_vntr_length]

    # Hexamer
    hexamer = 'CCCTCT'
    sva_hexamer_length = int(row['SVA_Hexamer'])  # Desired length of hexamer
    # Repeat the hexamer until we exceed or match the required length, then slice to get the exact length
    hexamer_seq = hexamer * (sva_hexamer_length // 6) + hexamer[:sva_hexamer_length % 6]

    # Transduction
    # Fetch a sequence using pysam
    with pysam.FastaFile(reference_fasta) as fasta_file:
        transduction = fasta_file.fetch(row['SRC_ref'], row['TD_beg'], row['TD_end'])

    # Final sequence
    seq = transduction + hexamer_seq + Alu_like_seq + vntr_sequence + SINE_R_seq + polyA1 + TSD

    return seq

def SVA__Hexamer_Alu_like_VNTR_SINE_R_POLYA_TD_POLYA_insertions(row, dict_consensus, reference_fasta, SVA_VNTR_list):
    # Get the L1_Seq from the dictionary
    SINE_R_seq = dict_consensus.get('SVA_SINE-R_Seq', '')
    Alu_like_seq = dict_consensus.get('SVA_Alu-like_Seq', '')

    # Poly A 1
    if pd.isna(row['PolyA_Length_1']) or row['PolyA_Length_1'] == 'NA':
        polyA_length1 = 0
    else:
        polyA_length1 = int(row['PolyA_Length_1'])
    polyA1 = 'A' * polyA_length1

    # Poly A 2
    if pd.isna(row['PolyA_Length_2']) or row['PolyA_Length_2'] == 'NA':
        polyA_length2 = 0 
    else:
        polyA_length2 = int(row['PolyA_Length_2'])
    polyA2 = 'A' * polyA_length2

    # TSD
    beg = int(row['beg'])
    tsd_len = int(row['TSD_Length'])
    start_tsd = beg - tsd_len
    # Fetch a sequence using pysam
    with pysam.FastaFile(reference_fasta) as fasta_file:
        TSD = fasta_file.fetch(row['#ref'], start_tsd, beg)

    # VNTR
    random_row = random.choice(SVA_VNTR_list).strip()  # Select random VNTR motif
    sva_vntr_length = int(row['SVA_VNTR_Length'])  # Desired length of VNTR
    
    # Generate VNTR sequence taking desired bases
    # In case VNTR is smaller than desired bases, repeat and wrap around the selected sequence
    vntr_sequence = ''
    while len(vntr_sequence) < sva_vntr_length:
        vntr_sequence += random_row
    vntr_sequence = vntr_sequence[:sva_vntr_length]

    # Hexamer
    hexamer = 'CCCTCT'
    sva_hexamer_length = int(row['SVA_Hexamer'])  # Desired length of hexamer
    # Repeat the hexamer until we exceed or match the required length, then slice to get the exact length
    hexamer_seq = hexamer * (sva_hexamer_length // 6) + hexamer[:sva_hexamer_length % 6]

    # Transduction
    # Fetch a sequence using pysam
    with pysam.FastaFile(reference_fasta) as fasta_file:
        transduction = fasta_file.fetch(row['SRC_ref'], row['TD_beg'], row['TD_end'])

    # Final sequence
    seq = hexamer_seq + Alu_like_seq + vntr_sequence + SINE_R_seq + polyA1 + transduction + polyA2 + TSD

    return seq

def generate_insertion_seq(row, motifs_file, reference_fasta, dict_consensus, SVA_VNTR_list):
    '''
    Global function that for each row of the df generates the determied sequence
    '''

    # VNTR
    if row['name'] == 'VNTR':
        return VNTR_insertions(row, motifs_file) 
    
    # DUPLICATIONS
    elif row['name'] == 'DUP':
        return DUP_insertions(row, reference_fasta)   
    
    # NUMT
    elif row['name'] == 'NUMT':
        return NUMT_insertions(row, dict_consensus)

    # INVERSE DUPLICATIONS  
    elif row['name'] == 'INV_DUP':
        seq = DUP_insertions(row, reference_fasta) 
        return inverse_sequence(seq)

    # ORPHAN
    elif row['name'] == 'orphan':
        return orphan_insertions(row, reference_fasta)   

    # Alu__FOR+POLYA
    elif row['name'] == 'Alu__FOR+POLYA':
        return Alu__FOR_POLYA_insertions(row, dict_consensus, reference_fasta)   

    # Alu__TRUN+FOR+POLYA
    elif row['name'] == 'Alu__TRUN+FOR+POLYA':
        return Alu__FOR_POLYA_insertions(row, dict_consensus, reference_fasta)   

    # L1__FOR+POLYA
    elif row['name'] == 'L1__FOR+POLYA':
        return L1__FOR_POLYA_insertions(row, dict_consensus, reference_fasta)   
    
    # L1__TRUN+FOR+POLYA
    elif row['name'] == 'L1__TRUN+FOR+POLYA':
        return L1__FOR_POLYA_insertions(row, dict_consensus, reference_fasta)   
    
    # L1__TD+FOR+POLYA
    elif row['name'] == 'L1__TD+FOR+POLYA':
        return L1__TD_FOR_POLYA_insertions(row, dict_consensus, reference_fasta)   

    # L1__FOR+POLYA+TD+POLYA
    elif row['name'] == 'L1__FOR+POLYA+TD+POLYA':
        return L1__FOR_POLYA_TD_POLYA_insertions(row, dict_consensus, reference_fasta)   

    # L1__TRUN+FOR+POLYA+TD+POLYA
    elif row['name'] == 'L1__TRUN+FOR+POLYA+TD+POLYA':
        return L1__FOR_POLYA_TD_POLYA_insertions(row, dict_consensus, reference_fasta)   

    # L1__TRUN+REV+BLUNT+FOR+POLYA
    elif row['name'] == 'L1__TRUN+REV+BLUNT+FOR+POLYA':
        return L1__TRUN_REV_BLUNT_FOR_POLYA_insertions(row, dict_consensus, reference_fasta) 

    # L1__TRUN+REV+BLUNT+FOR+POLYA+TD+POLYA
    elif row['name'] == 'L1__TRUN+REV+BLUNT+FOR+POLYA+TD+POLYA':
        return L1__TRUN_REV_BLUNT_FOR_POLYA_TD_POLYA_insertions(row, dict_consensus, reference_fasta) 
    
    # L1__TRUN+REV+DUP+FOR+POLYA
    elif row['name'] == 'L1__TRUN+REV+DUP+FOR+POLYA':
        return L1__TRUN_REV_DUP_FOR_POLYA_insertions(row, dict_consensus, reference_fasta) 

    # L1__TRUN+REV+DUP+FOR+POLYA+TD+POLYA
    elif row['name'] == 'L1__TRUN+REV+DUP+FOR+POLYA+TD+POLYA':
        return L1__TRUN_REV_DUP_FOR_POLYA_TD_POLYA_insertions(row, dict_consensus, reference_fasta) 

    # L1__TRUN+REV+DEL+FOR+POLYA
    elif row['name'] == 'L1__TRUN+REV+DEL+FOR+POLYA':
        return L1__TRUN_REV_DEL_FOR_POLYA_insertions(row, dict_consensus, reference_fasta) 
    
    # L1__REV+DEL+FOR+POLYA
    elif row['name'] == 'L1__REV+DEL+FOR+POLYA':
        return L1__TRUN_REV_DEL_FOR_POLYA_insertions(row, dict_consensus, reference_fasta) 

    # L1__TRUN+REV+DEL+FOR+POLYA+TD+POLYA
    elif row['name'] == 'L1__TRUN+REV+DEL+FOR+POLYA+TD+POLYA':
        return L1__TRUN_REV_DEL_FOR_POLYA_TD_POLYA_insertions(row, dict_consensus, reference_fasta) 

    # SVA__SINE-R+POLYA
    elif row['name'] == 'SVA__SINE-R+POLYA':
        return SVA__SINE_R_POLYA_insertions(row, dict_consensus, reference_fasta) 

    # SVA__VNTR+SINE-R+POLYA
    elif row['name'] == 'SVA__VNTR+SINE-R+POLYA':
        return SVA__VNTR_SINE_R_POLYA_insertions(row, dict_consensus, reference_fasta, SVA_VNTR_list)    

    # SVA__Alu-like+VNTR+SINE-R+POLYA
    elif row['name'] == 'SVA__Alu-like+VNTR+SINE-R+POLYA':
        return SVA__Alu_like_VNTR_SINE_R_POLYA_insertions(row, dict_consensus, reference_fasta, SVA_VNTR_list)    

    # SVA__MAST2+VNTR+SINE-R+POLYA
    elif row['name'] == 'SVA__MAST2+VNTR+SINE-R+POLYA':
        return SVA__MAST2_VNTR_SINE_R_POLYA_insertions(row, dict_consensus, reference_fasta, SVA_VNTR_list)  

    # SVA__TD+MAST2+VNTR+SINE-R+POLYA
    elif row['name'] == 'SVA__TD+MAST2+VNTR+SINE-R+POLYA':
        return SVA__TD_MAST2_VNTR_SINE_R_POLYA_insertions(row, dict_consensus, reference_fasta, SVA_VNTR_list)    

    # SVA__SINE-R+POLYA+TD+POLYA
    elif row['name'] == 'SVA__SINE-R+POLYA+TD+POLYA':
        return SVA__SINE_R_POLYA_TD_POLYA_insertions(row, dict_consensus, reference_fasta)   

    # SVA__VNTR+SINE-R+POLYA+TD+POLYA
    elif row['name'] == 'SVA__VNTR+SINE-R+POLYA+TD+POLYA':
        return SVA__VNTR_SINE_R_POLYA_TD_POLYA_insertions(row, dict_consensus, reference_fasta, SVA_VNTR_list)

    # SVA__Alu-like+VNTR+SINE-R+POLYA+TD+POLYA
    elif row['name'] == 'SVA__Alu-like+VNTR+SINE-R+POLYA+TD+POLYA':
        return SVA__Alu_like_VNTR_SINE_R_POLYA_TD_POLYA_insertions(row, dict_consensus, reference_fasta, SVA_VNTR_list)

    # SVA__MAST2+VNTR+SINE-R+POLYA+TD+POLYA
    elif row['name'] == 'SVA__MAST2+VNTR+SINE-R+POLYA+TD+POLYA':
        return SVA__MAST2_VNTR_SINE_R_POLYA_TD_POLYA_insertions(row, dict_consensus, reference_fasta, SVA_VNTR_list)  

    # SVA__Hexamer+Alu-like+VNTR+SINE-R+POLYA
    elif row['name'] == 'SVA__Hexamer+Alu-like+VNTR+SINE-R+POLYA':
        return SVA__Hexamer_Alu_like_VNTR_SINE_R_POLYA_insertions(row, dict_consensus, reference_fasta, SVA_VNTR_list)    

    # SVA__TD+Hexamer+Alu-like+VNTR+SINE-R+POLYA
    elif row['name'] == 'SVA__TD+Hexamer+Alu-like+VNTR+SINE-R+POLYA':
        return SVA__TD_Hexamer_Alu_like_VNTR_SINE_R_POLYA_insertions(row, dict_consensus, reference_fasta, SVA_VNTR_list)    
    
    # SVA__Hexamer+Alu-like+VNTR+SINE-R+POLYA+TD+POLYA
    elif row['name'] == 'SVA__Hexamer+Alu-like+VNTR+SINE-R+POLYA+TD+POLYA':
        return SVA__Hexamer_Alu_like_VNTR_SINE_R_POLYA_TD_POLYA_insertions(row, dict_consensus, reference_fasta, SVA_VNTR_list)   

    else:
        return 0
    
def reverse_complementary(seq):
    ''' 
    Function to create reverse complementary of a given sequence
    '''
    # initialize empty comp list to store complementary sequence
    comp = []
    seq = seq.upper()
    for base in seq:
        if base == 'A':
            comp.append('T')
        elif base == 'T':
            comp.append('A')
        elif base == 'G':
            comp.append('C')
        elif base == 'C':
            comp.append('G')
    # start counter with the length of the complementary sequence 
    counter = len(comp) - 1
    # empty list to start the reversee
    reverse = []
    # while the counter is not 0
    while counter >= 0:
        # take the last position of the complementary sequence
        base = comp[counter]
        # add it to the rerverse list
        reverse.append(base) 
        # jump to the previous position
        counter -= 1
    # join the list
    reverse_sequence = ''.join(reverse)   
    return reverse_sequence

def RC_insertion(df_insertions):
    '''
    Function to create the reverse complementary of those new generated sequences in the negative strand
    '''
    df_insertions['Sequence_Insertion'] = df_insertions.apply(
        lambda row: reverse_complementary(row['Sequence_Insertion']) if row['Strand'] == '-' else row['Sequence_Insertion'],
        axis=1
    )
    return df_insertions

def update_sequences(df_insertions):
    # Add a new column with 'Insertion' in each row
    df_insertions['Event_Type'] = 'Insertion'
    # Write all sequences in upper case
    df_insertions['Sequence_Insertion'] = df_insertions['Sequence_Insertion'].str.upper()

    # For those in - strand apply reverse complementary
    df_insertions = RC_insertion(df_insertions)

    return df_insertions

def main(consensus_path,probabilities_numbers_path,insertion_features_path,genome_wide_path,source_L1_path,source_SVA_path,motifs_path,SVA_VNTR_path,reference_fasta_path,num_events):
    print(f'File with consensus sequences: {consensus_path}')
    print(f'File with probabilities or number of events: {probabilities_numbers_path}')
    print(f'File with insertions features: {insertion_features_path}')
    print(f'File with genome-wide distribution: {genome_wide_path}')
    print(f'File with source genes for L1 transductions: {source_L1_path}')
    print(f'File with source genes for SVA transductions: {source_SVA_path}')
    print(f'File with VNTR motifs: {motifs_path}')
    print(f'File with reference genome: {reference_fasta_path}')
    print(f'File with SVA VNTR motifs: {SVA_VNTR_path}')
    print(f'Number of generated events: {num_events}')
    
    # Get consensus sequences
    consensus_dict = consensus_seqs(consensus_path)
    # Open SVAs VNTR motifs file
    SVA_VNTR_motifs = read_file_and_store_lines(SVA_VNTR_path)
    # Get number or probabilities of events
    df_insertions1 = probabilities_total_number(probabilities_numbers_path,num_events)
    # Process insetion features and generate dict with random numbers based on distributions for each feature of every event
    dict_random = process_insertion_features_random_numbers(insertion_features_path,num_events,consensus_dict)
    # Add chromosome and start position
    df_insertions2 = add_beg_end_columns(df_insertions1,genome_wide_path)
    # Add the features of each event
    df_insertions3 = add_elements_columns(dict_random,df_insertions2)
    # Add SVA events additional details
    df_insertions4 = add_SVA_info(df_insertions3, consensus_dict)
    # Update df
    df_insertions5 = update_dataframe(df_insertions4,consensus_dict)
    # Add souce gene information for transduction events
    df_insertions6 = add_source_gene_info(df_insertions5,source_L1_path,source_SVA_path)
    # Generate the insertion sequence
    df_insertions6['Sequence_Insertion'] = df_insertions6.apply(lambda row: generate_insertion_seq(row, motifs_path, reference_fasta_path, consensus_dict, SVA_VNTR_motifs), axis=1)
    # Update the insertion sequence
    df_insertions6 = update_sequences(df_insertions6)
    # Save the output
    df_insertions6.to_csv('Insertions_table.tsv', sep='\t', index=False)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Generate insertion sequences returned in tsv file.')
    parser.add_argument('consensus_path', type=str, help='Path to file with consensus sequences.')
    parser.add_argument('probabilities_numbers_path', type=str, help='Path to the TSV file probabilities or defined number of events.')
    parser.add_argument('insertion_features_path', type=str, help='Path to the TSV file containing insertion features of events.')
    parser.add_argument('genome_wide_path', type=str, help='Path to the TSV file containing genome-wide distribution of events.')
    parser.add_argument('source_L1_path', type=str, help='Path to the TSV file containing loci for LINE-1 transductions.')
    parser.add_argument('source_SVA_path', type=str, help='Path to the TSV file containing loci for SVA transductions.')
    parser.add_argument('motifs_path', type=str, help='Path to the txt file containing VNTR motifs and their proportion.')
    parser.add_argument('SVA_VNTR_path', type=str, help='Path to the txt file containing SVA VNTR motifs.')
    parser.add_argument('reference_fasta_path', type=str, help='Path to file with reference genome.')
    parser.add_argument('--num_events', type=int, default=100, help='Number of events to sample (optional, just in case of providing probabilities).')
    args = parser.parse_args()
    main(args.consensus_path, args.probabilities_numbers_path, args.insertion_features_path, args.genome_wide_path, args.source_L1_path, args.source_SVA_path, args.motifs_path, args.SVA_VNTR_path, args.reference_fasta_path, args.num_events)
