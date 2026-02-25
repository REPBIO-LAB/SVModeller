# SVModeller - Module 2

# Model probabilities, insertion features, and genome wide distribution to generate new events

# Input:
# - Genome-Wide Distribution 
# - Insertion Features 
# - Event Probabilities or Number of each event to simulate 
# - OPTIONAL: Just in case of providing probabilities, total number of events to simulate
# - Table with source loci to LINE-1 transductions
# - Table with source loci to SVA transductions
# - Consensus sequences
# - Reference genome 
# - List of VNTR motifs
# - List of SVA VNTR motifs 

# Output:
# - New insertion events sequences with their corresponding features (Insertions_table.tsv)
# - OPTIONAL: Variant Calling File (VCF) with insertion data

# Developers
# SVModeller has been developed by Ismael Vera-Munoz (orcid.org/0009-0009-2860-378X) at the Repetitive DNA Biology (REPBIO) Lab at the Centre for Genomic Regulation (CRG) (Barcelona 2024-2026)

# License
# SVModeller is distributed under the AGPL-3.0.

import argparse
import numpy as np
import pandas as pd
import distfit
from distfit import distfit
import random
import pysam
from GAPI import formats
import os
import datetime
import warnings

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
        dist.fit_transform(value, verbose = 0)
        random_numbers = dist.generate(num_samples)
        random_numbers = random_numbers.astype(int)
        random_numbers_dict[key] = random_numbers
    return random_numbers_dict

def remove_negative_values(values):
    ''' 
    Function to remove negative values from the list
    '''
    return [value for value in values if value > 0]

def filter_sd(dict_mutations, key_list):
    '''
    Remove values outside +-2 standard deviations from the mean
    Only applied to numeric values
    '''
    for key in key_list:
        values = dict_mutations.get(key, [])

        # Keep only numeric values
        numeric_values = [v for v in values if isinstance(v, (int, float))]

        # Skip if not enough numeric values
        if len(numeric_values) < 2:
            continue

        mean = np.mean(numeric_values)
        std_dev = np.std(numeric_values, ddof=1)

        upper_lim = mean + 2 * std_dev
        lower_lim = mean - 2 * std_dev

        # Filter only numeric entries, keep non-numeric untouched
        filtered = []
        for v in values:
            if isinstance(v, (int, float)):
                if lower_lim <= v <= upper_lim:
                    filtered.append(v)
            else:
                filtered.append(v)

        dict_mutations[key] = filtered

    return dict_mutations

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
        array = np.array(dict_insertion_features[key])
        array_dict[key] = array

    number_random_events = num_events * 3
    random_numbers_dict = generate_random_numbers(array_dict, number_random_events)

    # Remove negative values
    for key, values in random_numbers_dict.items():
        random_numbers_dict[key] = remove_negative_values(values)

    # Filter values outside ±2 standard deviations
    keys_filter = list(random_numbers_dict.keys())
    random_numbers_dict = filter_sd(random_numbers_dict, keys_filter)

    # Filter largest DEL values
    random_numbers_dict = filter_DEL(random_numbers_dict)

    # Filter FOR values larger than consensus sequences
    random_numbers_dict = filter_FOR(random_numbers_dict, dict_consensus)

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
    Function to generate VNTR sequences by selecting a random row from the provided file
    and using its 'Complete_Sequence', 'Start', 'VNTR_Num_Motifs', and 'VNTR_Motifs'.
    '''

    # Step 1: Read the motifs file
    motifs_df = pd.read_csv(motifs_file, sep='\t')  # Reads the .tsv file containing motif information
    
    # Step 2: Randomly select a row from the motifs file
    random_row = motifs_df.sample(n=1).iloc[0]
    
    # Step 3: Get values from the randomly selected row
    complete_sequence = random_row['Complete_Sequence']
    start_position = random_row['Start']
    vntr_num_motifs = int(random_row['VNTR_Num_Motifs'])
    vntr_motifs = random_row['VNTR_Motifs'].split(',')
    # Step 4: Update the row with the selected values
    row['Sequence_Insertion'] = complete_sequence
    row['Length'] = len(complete_sequence)  # Set the Length as the length of the Complete_Sequence
    row['VNTR_Num_Motifs'] = vntr_num_motifs  # Update the VNTR_Num_Motifs
    row['VNTR_Motifs'] = vntr_motifs  # Update the VNTR_Motifs list
    row['Start'] = start_position  # Update the Start position
    
    # Step 5: Return the Complete_Sequence as the sequence and the VNTR_Motifs
    sequence = complete_sequence  # Now the sequence is just the Complete_Sequence
    return sequence, vntr_motifs, vntr_num_motifs

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
        rev_seq = reverse_complementary(rev_seq)  # Apply reverse complementary to REV sequence
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
        rev_seq = reverse_complementary(rev_seq)  # Apply reverse complementary to REV sequence
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
        rev_seq = reverse_complementary(rev_seq)  # Apply reverse complementary to REV sequence
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
        rev_seq = reverse_complementary(rev_seq)  # Apply reverse complementary to REV sequence
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
        rev_seq = reverse_complementary(rev_seq)  # Apply reverse complementary to REV sequence
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
        rev_seq = reverse_complementary(rev_seq)  # Apply reverse complementary to REV sequence
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
        sequence, selected_motifs, vntr_num_motifs = VNTR_insertions(row, motifs_file)
        return sequence, selected_motifs, vntr_num_motifs

    # DUPLICATIONS
    elif row['name'] == 'DUP':
        return DUP_insertions(row, reference_fasta), 0, 0
    
    # NUMT
    elif row['name'] == 'NUMT':
        return NUMT_insertions(row, dict_consensus), 0, 0

    # INVERSE DUPLICATIONS  
    elif row['name'] == 'INV_DUP':
        seq = DUP_insertions(row, reference_fasta)
        return reverse_complementary(seq), 0, 0

    # ORPHAN
    elif row['name'] == 'orphan':
        return orphan_insertions(row, reference_fasta)   , 0, 0

    # Alu__FOR+POLYA
    elif row['name'] == 'Alu__FOR+POLYA':
        return Alu__FOR_POLYA_insertions(row, dict_consensus, reference_fasta)   , 0, 0

    # Alu__TRUN+FOR+POLYA
    elif row['name'] == 'Alu__TRUN+FOR+POLYA':
        return Alu__FOR_POLYA_insertions(row, dict_consensus, reference_fasta)   , 0, 0

    # L1__FOR+POLYA
    elif row['name'] == 'L1__FOR+POLYA':
        return L1__FOR_POLYA_insertions(row, dict_consensus, reference_fasta)   , 0, 0
    
    # L1__TRUN+FOR+POLYA
    elif row['name'] == 'L1__TRUN+FOR+POLYA':
        return L1__FOR_POLYA_insertions(row, dict_consensus, reference_fasta) , 0  , 0
    
    # L1__TD+FOR+POLYA
    elif row['name'] == 'L1__TD+FOR+POLYA':
        return L1__TD_FOR_POLYA_insertions(row, dict_consensus, reference_fasta)   , 0, 0

    # L1__FOR+POLYA+TD+POLYA
    elif row['name'] == 'L1__FOR+POLYA+TD+POLYA':
        return L1__FOR_POLYA_TD_POLYA_insertions(row, dict_consensus, reference_fasta)   , 0, 0

    # L1__TRUN+FOR+POLYA+TD+POLYA
    elif row['name'] == 'L1__TRUN+FOR+POLYA+TD+POLYA':
        return L1__FOR_POLYA_TD_POLYA_insertions(row, dict_consensus, reference_fasta)   , 0, 0

    # L1__TRUN+REV+BLUNT+FOR+POLYA
    elif row['name'] == 'L1__TRUN+REV+BLUNT+FOR+POLYA':
        return L1__TRUN_REV_BLUNT_FOR_POLYA_insertions(row, dict_consensus, reference_fasta) , 0, 0

    # L1__TRUN+REV+BLUNT+FOR+POLYA+TD+POLYA
    elif row['name'] == 'L1__TRUN+REV+BLUNT+FOR+POLYA+TD+POLYA':
        return L1__TRUN_REV_BLUNT_FOR_POLYA_TD_POLYA_insertions(row, dict_consensus, reference_fasta) , 0, 0
    
    # L1__TRUN+REV+DUP+FOR+POLYA
    elif row['name'] == 'L1__TRUN+REV+DUP+FOR+POLYA':
        return L1__TRUN_REV_DUP_FOR_POLYA_insertions(row, dict_consensus, reference_fasta) , 0, 0

    # L1__TRUN+REV+DUP+FOR+POLYA+TD+POLYA
    elif row['name'] == 'L1__TRUN+REV+DUP+FOR+POLYA+TD+POLYA':
        return L1__TRUN_REV_DUP_FOR_POLYA_TD_POLYA_insertions(row, dict_consensus, reference_fasta) , 0, 0

    # L1__TRUN+REV+DEL+FOR+POLYA
    elif row['name'] == 'L1__TRUN+REV+DEL+FOR+POLYA':
        return L1__TRUN_REV_DEL_FOR_POLYA_insertions(row, dict_consensus, reference_fasta) , 0, 0
    
    # L1__REV+DEL+FOR+POLYA
    elif row['name'] == 'L1__REV+DEL+FOR+POLYA':
        return L1__TRUN_REV_DEL_FOR_POLYA_insertions(row, dict_consensus, reference_fasta) , 0, 0

    # L1__TRUN+REV+DEL+FOR+POLYA+TD+POLYA
    elif row['name'] == 'L1__TRUN+REV+DEL+FOR+POLYA+TD+POLYA':
        return L1__TRUN_REV_DEL_FOR_POLYA_TD_POLYA_insertions(row, dict_consensus, reference_fasta) , 0, 0

    # SVA__SINE-R+POLYA
    elif row['name'] == 'SVA__SINE-R+POLYA':
        return SVA__SINE_R_POLYA_insertions(row, dict_consensus, reference_fasta) , 0, 0

    # SVA__VNTR+SINE-R+POLYA
    elif row['name'] == 'SVA__VNTR+SINE-R+POLYA':
        return SVA__VNTR_SINE_R_POLYA_insertions(row, dict_consensus, reference_fasta, SVA_VNTR_list)   , 0 , 0

    # SVA__Alu-like+VNTR+SINE-R+POLYA
    elif row['name'] == 'SVA__Alu-like+VNTR+SINE-R+POLYA':
        return SVA__Alu_like_VNTR_SINE_R_POLYA_insertions(row, dict_consensus, reference_fasta, SVA_VNTR_list)    , 0, 0

    # SVA__MAST2+VNTR+SINE-R+POLYA
    elif row['name'] == 'SVA__MAST2+VNTR+SINE-R+POLYA':
        return SVA__MAST2_VNTR_SINE_R_POLYA_insertions(row, dict_consensus, reference_fasta, SVA_VNTR_list)  , 0, 0

    # SVA__TD+MAST2+VNTR+SINE-R+POLYA
    elif row['name'] == 'SVA__TD+MAST2+VNTR+SINE-R+POLYA':
        return SVA__TD_MAST2_VNTR_SINE_R_POLYA_insertions(row, dict_consensus, reference_fasta, SVA_VNTR_list)  , 0  , 0

    # SVA__SINE-R+POLYA+TD+POLYA
    elif row['name'] == 'SVA__SINE-R+POLYA+TD+POLYA':
        return SVA__SINE_R_POLYA_TD_POLYA_insertions(row, dict_consensus, reference_fasta)   , 0, 0

    # SVA__VNTR+SINE-R+POLYA+TD+POLYA
    elif row['name'] == 'SVA__VNTR+SINE-R+POLYA+TD+POLYA':
        return SVA__VNTR_SINE_R_POLYA_TD_POLYA_insertions(row, dict_consensus, reference_fasta, SVA_VNTR_list), 0, 0

    # SVA__Alu-like+VNTR+SINE-R+POLYA+TD+POLYA
    elif row['name'] == 'SVA__Alu-like+VNTR+SINE-R+POLYA+TD+POLYA':
        return SVA__Alu_like_VNTR_SINE_R_POLYA_TD_POLYA_insertions(row, dict_consensus, reference_fasta, SVA_VNTR_list), 0, 0

    # SVA__MAST2+VNTR+SINE-R+POLYA+TD+POLYA
    elif row['name'] == 'SVA__MAST2+VNTR+SINE-R+POLYA+TD+POLYA':
        return SVA__MAST2_VNTR_SINE_R_POLYA_TD_POLYA_insertions(row, dict_consensus, reference_fasta, SVA_VNTR_list)  , 0, 0

    # SVA__Hexamer+Alu-like+VNTR+SINE-R+POLYA
    elif row['name'] == 'SVA__Hexamer+Alu-like+VNTR+SINE-R+POLYA':
        return SVA__Hexamer_Alu_like_VNTR_SINE_R_POLYA_insertions(row, dict_consensus, reference_fasta, SVA_VNTR_list)    , 0, 0

    # SVA__TD+Hexamer+Alu-like+VNTR+SINE-R+POLYA
    elif row['name'] == 'SVA__TD+Hexamer+Alu-like+VNTR+SINE-R+POLYA':
        return SVA__TD_Hexamer_Alu_like_VNTR_SINE_R_POLYA_insertions(row, dict_consensus, reference_fasta, SVA_VNTR_list)    , 0, 0
    
    # SVA__Hexamer+Alu-like+VNTR+SINE-R+POLYA+TD+POLYA
    elif row['name'] == 'SVA__Hexamer+Alu-like+VNTR+SINE-R+POLYA+TD+POLYA':
        return SVA__Hexamer_Alu_like_VNTR_SINE_R_POLYA_TD_POLYA_insertions(row, dict_consensus, reference_fasta, SVA_VNTR_list)   , 0 , 0

    else:
        return '', 0, 0 

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

# Function to process each row's 'name' value
def process_name(name):
    if '__' in name:  # Case when there is '__' in the name
        fam_n, conformation = name.split('__', 1)
        itype_n = 'partnered' if 'TD' in conformation else 'solo'
    else:  # Case when there is no '__' in the name
        fam_n, conformation = np.nan, np.nan
        itype_n = name
    
    return pd.Series([itype_n, conformation, fam_n])

# Function to generate the HEXAMER_SEQ column
def generate_hexamer_seq(row):
    # Check if 'FAM_N' contains 'SVA' and 'CONFORMATION' contains 'Hexamer'
    if pd.notna(row['FAM_N']) and 'SVA' in row['FAM_N'] and pd.notna(row['CONFORMATION']) and 'Hexamer' in row['CONFORMATION']:
        hexamer = 'CCCTCT'
        try:
            # Assuming HEXAMER_LEN is a column with the desired hexamer length
            sva_hexamer_length = int(row['HEXAMER_LEN'])  # Desired length of hexamer
            # Generate the hexamer sequence
            hexamer_seq = hexamer * (sva_hexamer_length // 6) + hexamer[:sva_hexamer_length % 6]
            return hexamer_seq
        except (ValueError, TypeError):
            return np.nan  # Return NaN if there's an issue with the length
    return np.nan  # Return NaN if conditions are not met

# Function to generate TSD sequence from the reference fasta file
def generate_tsd_seq(row, reference_fasta):
    # Check if 'beg' and 'TSD_LEN' are valid
    if pd.notna(row['beg']) and pd.notna(row['TSD_LEN']):
        try:
            beg = int(row['beg'])
            tsd_len = int(row['TSD_LEN'])
            start_tsd = beg - tsd_len

            # Fetch the sequence from the reference FASTA file using pysam
            with pysam.FastaFile(reference_fasta) as fasta_file:
                TSD = fasta_file.fetch(row['#ref'], start_tsd, beg)  # Extract the sequence
                return TSD.upper()
        except (ValueError, TypeError, KeyError):
            return np.nan  # Return NaN if there's an issue with the length or values
    return np.nan  # Return NaN if either 'beg' or 'TSD_LEN' is missing or invalid

# Function to generate sequence for 3PRIME_TD
def generate_3prime_td_seq(row, reference_fasta):
    # Fetch the sequence for 3PRIME_TD_SEQ if there's a value in 3PRIME_NB_TD
    if pd.notna(row['3PRIME_NB_TD']):
        try:
            with pysam.FastaFile(reference_fasta) as fasta_file:
                transduction = fasta_file.fetch(row['SRC_ref'], row['TD_beg'], row['TD_end'])
                return transduction
        except (KeyError, ValueError):
            return np.nan
    return np.nan

# Function to generate sequence for 5PRIME_TD
def generate_5prime_td_seq(row, reference_fasta):
    # Fetch the sequence for 5PRIME_TD_SEQ if there's a value in 5PRIME_NB_TD
    if pd.notna(row['5PRIME_NB_TD']):
        try:
            with pysam.FastaFile(reference_fasta) as fasta_file:
                transduction = fasta_file.fetch(row['SRC_ref'], row['TD_beg'], row['TD_end'])
                return transduction
        except (KeyError, ValueError):
            return np.nan
    return np.nan

# Function to generate sequence for orphan TD
def generate_orphan_seq(row, reference_fasta):
    # Fetch the sequence for orphan TD if there's a value in ORPHAN_TD_COORD
    if pd.notna(row['ORPHAN_TD_COORD']):
        try:
            with pysam.FastaFile(reference_fasta) as fasta_file:
                transduction = fasta_file.fetch(row['SRC_ref'], row['TD_beg'], row['TD_end'])
                return transduction
        except (KeyError, ValueError):
            return np.nan
    return np.nan

# Function to generate POLYA_LEN column
def generate_polya_len(row):
    # Remove the decimal by converting the value to an integer before adding to the string
    poly_a_len_1 = str(int(row['PolyA_Length_1'])) if pd.notna(row['PolyA_Length_1']) else ''
    poly_a_len_2 = str(int(row['PolyA_Length_2'])) if pd.notna(row['PolyA_Length_2']) else ''
    
    if poly_a_len_1 and poly_a_len_2:
        return poly_a_len_1 + ',' + poly_a_len_2
    elif poly_a_len_1:
        return poly_a_len_1
    else:
        return np.nan  # Return NaN if both are missing

# Function to generate POLYA_SEQ column
def generate_polya_seq(row):
    strand = row['STRAND'] if pd.notna(row['STRAND']) else '+'
    
    # Determine the PolyA sequences
    seq_1 = 'A' * int(row['PolyA_Length_1']) if pd.notna(row['PolyA_Length_1']) else ''
    seq_2 = 'A' * int(row['PolyA_Length_2']) if pd.notna(row['PolyA_Length_2']) else ''
    
    if strand == '-':
        seq_1 = seq_1.replace('A', 'T') if seq_1 else ''
        seq_2 = seq_2.replace('A', 'T') if seq_2 else ''
    
    # Join the sequences with a comma if both exist
    if seq_1 and seq_2:
        return seq_1 + ',' + seq_2
    elif seq_1:
        return seq_1
    elif seq_2:
        return seq_2
    return np.nan  # Return NaN if no sequences are available

def process_vntr_motifs(df):
    # Iterate over each element in the 'Selected_VNTR_Motifs' column
    for i in range(len(df)):
        #motif = df.at[i, 'Selected_VNTR_Motifs']
        motif = df.iloc[i]['Selected_VNTR_Motifs']
        if motif == '0':
            df.at[i, 'Selected_VNTR_Motifs'] = np.nan  # Replace '0' with NaN
        elif isinstance(motif, str):
            # Remove '[' and ']' from the string
            df.at[i, 'Selected_VNTR_Motifs'] = motif.replace("[", "").replace("'[", "").replace("]", "").replace("]'", "").strip()
    
    return df

# Main function
def df_VCF(df, reference_fasta):
    # Create a copy of the DataFrame
    df_copy = df.copy()
    df_copy = process_vntr_motifs(df_copy)
    df_copy['Selected_VNTR_Motifs'] = df_copy['Selected_VNTR_Motifs'].str.replace("'", "")
    # Rename columns based on your requirements
    df_copy = df_copy.rename(columns={
        'Length': 'INS_LEN',
        'VNTR_Num_Motifs': 'NB_MOTIFS',
        'SVA_Hexamer': 'HEXAMER_LEN',
        'TSD_Length': 'TSD_LEN',
        'TD_5': '5PRIME_TD_LEN',
        'TD_3': '3PRIME_TD_LEN',
        'Strand': 'STRAND',
        'TD_orphan_Length': 'ORPHAN_TD_LEN',
        'Selected_VNTR_Motifs': 'MOTIFS'
    })

    # Convert columns to integers to remove decimals (i.e., .0)
    df_copy['NB_MOTIFS'] = df_copy['NB_MOTIFS'].apply(lambda x: int(x) if pd.notna(x) else np.nan).astype('Int64')  # Use 'Int64' to keep NaNs
    df_copy['HEXAMER_LEN'] = df_copy['HEXAMER_LEN'].apply(lambda x: int(x) if pd.notna(x) else np.nan).astype('Int64')  # Use 'Int64' to keep NaNs
    df_copy['TSD_LEN'] = df_copy['TSD_LEN'].apply(lambda x: int(x) if pd.notna(x) else np.nan).astype('Int64')  # Use 'Int64' to keep NaNs
    df_copy['5PRIME_TD_LEN'] = df_copy['5PRIME_TD_LEN'].apply(lambda x: int(x) if pd.notna(x) else np.nan).astype('Int64')  # Use 'Int64' to keep NaNs
    df_copy['3PRIME_TD_LEN'] = df_copy['3PRIME_TD_LEN'].apply(lambda x: int(x) if pd.notna(x) else np.nan).astype('Int64')  # Use 'Int64' to keep NaNs
    
    # Add the 'ID' column with the format SV1, SV2, SV3...
    df_copy['ID'] = ['INS_' + str(i + 1) for i in range(len(df_copy))]

    # Apply the processing function to the 'name' column and split it into three new columns
    df_copy[['ITYPE_N', 'CONFORMATION', 'FAM_N']] = df_copy['name'].apply(process_name)
    
    # For rows where there is no 'CONFORMATION' (e.g., VNTR), set 'CONFORMATION' to NaN
    df_copy['CONFORMATION'] = df_copy['CONFORMATION'].replace('', np.nan)

    # Remove the 'name' column at the end
    df_copy = df_copy.drop(columns=['name'])
    
    # Apply the hexamer sequence generation logic
    df_copy['HEXAMER_SEQ'] = df_copy.apply(generate_hexamer_seq, axis=1)
    
    # Apply the TSD sequence generation logic
    df_copy['TSD_SEQ'] = df_copy.apply(lambda row: generate_tsd_seq(row, reference_fasta), axis=1)
    
    # Create the 3PRIME_NB_TD and 5PRIME_NB_TD columns based on the conditions
    df_copy['3PRIME_NB_TD'] = df_copy['3PRIME_TD_LEN'].apply(lambda x: 1 if pd.notna(x) else np.nan)
    df_copy['5PRIME_NB_TD'] = df_copy['5PRIME_TD_LEN'].apply(lambda x: 1 if pd.notna(x) else np.nan)
    
    # Create the 3PRIME_TD_COORD, 5PRIME_TD_COORD and ORPHAN_TD_COORD columns based on the conditions
    df_copy['3PRIME_TD_COORD'] = df_copy['SRC_identifier'].where(df_copy['3PRIME_NB_TD'].notna(), np.nan)
    df_copy['5PRIME_TD_COORD'] = df_copy['SRC_identifier'].where(df_copy['5PRIME_NB_TD'].notna(), np.nan)
    df_copy['ORPHAN_TD_COORD'] = df_copy['SRC_identifier'].where(df_copy['ITYPE_N'] == 'orphan')

    # Remove the SRC_identifier column
    df_copy = df_copy.drop(columns=['SRC_identifier'])

    # Apply the transduction sequence generation logic for 3PRIME_TD_SEQ and 5PRIME_TD_SEQ
    df_copy['3PRIME_TD_SEQ'] = df_copy.apply(lambda row: generate_3prime_td_seq(row, reference_fasta), axis=1)
    df_copy['5PRIME_TD_SEQ'] = df_copy.apply(lambda row: generate_5prime_td_seq(row, reference_fasta), axis=1)
    df_copy['ORPHAN_TD_SEQ'] = df_copy.apply(lambda row: generate_orphan_seq(row, reference_fasta), axis=1)  

    # Convert the 3PRIME_TD_SEQ and 5PRIME_TD_SEQ columns to uppercase
    df_copy['3PRIME_TD_SEQ'] = df_copy['3PRIME_TD_SEQ'].apply(lambda x: x.upper() if pd.notna(x) else np.nan)
    df_copy['5PRIME_TD_SEQ'] = df_copy['5PRIME_TD_SEQ'].apply(lambda x: x.upper() if pd.notna(x) else np.nan)
    df_copy['ORPHAN_TD_SEQ'] = df_copy['ORPHAN_TD_SEQ'].apply(lambda x: x.upper() if pd.notna(x) else np.nan)

    # Apply the POLYA_LEN generation logic
    df_copy['POLYA_LEN'] = df_copy.apply(generate_polya_len, axis=1)

    # Apply the POLYA_SEQ generation logic
    df_copy['POLYA_SEQ'] = df_copy.apply(generate_polya_seq, axis=1)

    return df_copy

def create_vcf_file(df, reference_fasta, chromosome_length):
    # Create the chr_length dictionary using formats
    chr_length = formats.chrom_lengths_index(chromosome_length)

    # Get date of creation
    current_date = datetime.datetime.now().strftime("%Y%m%d")
    reference_name = os.path.basename(reference_fasta)

    # Order chromosomes
    def sort_chromosomes(contig):
        contig_clean = contig.lower().replace("chr", "")
        special_order = {"x": 23, "y": 24, "m": 25, "mt": 25}

        if contig_clean.isdigit():
            return (0, int(contig_clean))

        if contig_clean in special_order:
            return (0, special_order[contig_clean])

        return (1, contig_clean)

    contigs = sorted(df['#ref'].unique(), key=sort_chromosomes)

    # Open a VCF file to write to
    with open('VCF_Insertions_SVModeller.vcf', 'w') as vcf_file:
        # Write VCF header
        vcf_file.write("##fileformat=VCFv4.2\n")
        vcf_file.write(f"##fileDate={current_date}\n")
        vcf_file.write("##source=SVModeller\n")
        vcf_file.write(f"##reference={reference_name}\n")

        # Adding contig length
        for contig in contigs:
            if contig in chr_length:
                contig_length = chr_length[contig]
                vcf_file.write(f"##contig=<ID={contig},assembly=None,length={contig_length},species=None>\n")
            else:
                vcf_file.write(f"##contig=<ID={contig},assembly=None,length=None,species=None>\n")

        vcf_file.write('##INFO=<ID=SVTYPE,Number=1,Type=String,Description="Type of structural variant">\n')
        vcf_file.write("##INFO=<ID=INS_LEN,Type=float,Description=Total length of the insertion>\n")
        vcf_file.write("##INFO=<ID=TSD_LEN,Type=float,Description=Length of the TSD>\n")
        vcf_file.write("##INFO=<ID=5PRIME_TD_LEN,Type=float,Description=Length of the 5' TD>\n")
        vcf_file.write("##INFO=<ID=3PRIME_TD_LEN,Type=float,Description=Length of the 3' TD>\n")
        vcf_file.write("##INFO=<ID=NB_MOTIFS,Type=float,Description=Number of VNTR motifs>\n")
        vcf_file.write("##INFO=<ID=MOTIFS,Type=str,Description=VNTR selected motifs>\n")
        vcf_file.write("##INFO=<ID=HEXAMER_LEN,Type=float,Description=Length of the SVA hexamer>\n")
        vcf_file.write("##INFO=<ID=SVA_VNTR_Length,Type=float,Description=Length of the SVA VNTR>\n")
        vcf_file.write("##INFO=<ID=ITYPE_N,Type=float,Description=Type of insertion>\n")
        vcf_file.write("##INFO=<ID=CONFORMATION,Type=float,Description=Conformation of the insertion>\n")
        vcf_file.write("##INFO=<ID=FAM_N,Type=float,Description=Family of the insertion>\n")
        vcf_file.write("##INFO=<ID=HEXAMER_SEQ,Type=float,Description=Sequence of SVA hexamer>\n")
        vcf_file.write("##INFO=<ID=TSD_SEQ,Type=float,Description=Length of the poly A tail (or first poly A in case of more than 1)>\n")
        vcf_file.write("##INFO=<ID=3PRIME_NB_TD,Type=float,Description=Number of 3' TD>\n")
        vcf_file.write("##INFO=<ID=5PRIME_NB_TD,Type=float,Description=Number of 5' TD>\n")
        vcf_file.write("##INFO=<ID=3PRIME_TD_COORD,Type=float,Description=Coordinates of 3' TD>\n")
        vcf_file.write("##INFO=<ID=5PRIME_TD_COORD,Type=float,Description=Coordinates of 5' TD>\n")
        vcf_file.write("##INFO=<ID=3PRIME_TD_SEQ,Type=float,Description=Sequence of 3' TD>\n")
        vcf_file.write("##INFO=<ID=5PRIME_TD_SEQ,Type=float,Description=Sequence of 5' TD>\n")
        vcf_file.write("##INFO=<ID=POLYA_LEN,Type=float,Description=Length of the poly A tail>\n")
        vcf_file.write("##INFO=<ID=POLYA_SEQ,Type=float,Description=Sequence of the poly A tail>\n")
        vcf_file.write("##INFO=<ID=FOR,Type=float,Description=Length of conformation forward part of the event>\n")
        vcf_file.write("##INFO=<ID=TRUN,Type=float,Description=Length of the conformation truncated part of the event>\n")
        vcf_file.write("##INFO=<ID=REV,Type=float,Description=Length of the conformation reverse part of the event>\n")
        vcf_file.write("##INFO=<ID=DEL,Type=float,Description=Length of the conformation deleted part of the event>\n")
        vcf_file.write("##INFO=<ID=DUP,Type=float,Description=Length of the conformation duplicated part of the event>\n")
        vcf_file.write("##INFO=<ID=STRAND,Type=str,Description=Strand (+ or -) of the event>\n")
        vcf_file.write("##INFO=<ID=ORPHAN_TD_LEN,Type=float,Description=Length of the orphan transduction>\n")
        vcf_file.write("##INFO=<ID=ORPHAN_TD_COORD,Type=str,Description=Coordinates of orphan transduction>\n")
        vcf_file.write("##INFO=<ID=ORPHAN_TD_SEQ,Type=str,Description=Sequence of orphan transduction>\n")

        vcf_file.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n")

        with pysam.FastaFile(reference_fasta) as fasta_file:

            for _, row in df.iterrows():
                chrom = row['#ref']
                pos = row['beg']
                event_id = row['ID']
                alt = row['Sequence_Insertion']
                qual = '.'
                filter_val = '.'

                beg = int(pos)

                # Fetch reference base
                ref = fasta_file.fetch(chrom, beg - 1, beg).upper()

                info_fields = []
                for col in ['ITYPE_N', 'INS_LEN', 'STRAND', 'FAM_N', 'CONFORMATION',
                            'FOR', 'TRUN', 'REV', 'DEL', 'DUP', 'TSD_LEN', 'TSD_SEQ',
                            '3PRIME_NB_TD', '5PRIME_NB_TD', '5PRIME_TD_LEN',
                            '3PRIME_TD_LEN', '5PRIME_TD_COORD', '3PRIME_TD_COORD',
                            '3PRIME_TD_SEQ', '5PRIME_TD_SEQ', 'NB_MOTIFS',
                            'MOTIFS', 'HEXAMER_LEN', 'SVA_VNTR_Length',
                            'HEXAMER_SEQ', 'ORPHAN_TD_LEN', 'ORPHAN_TD_COORD',
                            'ORPHAN_TD_SEQ', 'POLYA_LEN', 'POLYA_SEQ']:
                    
                    value = row[col]
                    if pd.notna(value):
                        info_fields.append(f"{col}={value}")
                
                info_fields.append("SVTYPE=INS")

                info = ";".join(info_fields)

                vcf_file.write(f"{chrom}\t{pos}\t{event_id}\t{ref}\t{alt}\t{qual}\t{filter_val}\t{info}\n")

    with open('VCF_Insertions_SVModeller.vcf', 'r') as file:
        content = file.read()

    content = content.replace("5PRIME_TD_LEN", "TD_LEN_5PRIME")
    content = content.replace("3PRIME_TD_LEN", "TD_LEN_3PRIME")
    content = content.replace("3PRIME_NB_TD", "NB_TD_3PRIME")
    content = content.replace("5PRIME_NB_TD", "NB_TD_5PRIME")
    content = content.replace("3PRIME_TD_COORD", "TD_COORD_3PRIME")
    content = content.replace("5PRIME_TD_COORD", "TD_COORD_5PRIME")
    content = content.replace("3PRIME_TD_SEQ", "TD_SEQ_3PRIME")
    content = content.replace("5PRIME_TD_SEQ", "TD_SEQ_5PRIME")
    content = content.replace("INS_LEN", "SVLEN")

    with open('VCF_Insertions_SVModeller.vcf', 'w') as file:
        file.write(content)

    print("VCF file created successfully.")

# Remove FutureWarnings
warnings.simplefilter(action='ignore', category=FutureWarning)
warnings.simplefilter(action='ignore', category=UserWarning)

def main(consensus_path, probabilities_numbers_path, insertion_features_path, genome_wide_path, source_L1_path, source_SVA_path, motifs_path, SVA_VNTR_path, reference_fasta_path, chromosome_length_path, num_events, apply_VCF):
    print(f'File with consensus sequences: {consensus_path}')
    print(f'File with probabilities or number of events: {probabilities_numbers_path}')
    print(f'File with insertions features: {insertion_features_path}')
    print(f'File with genome-wide distribution: {genome_wide_path}')
    print(f'File with source genes for L1 transductions: {source_L1_path}')
    print(f'File with source genes for SVA transductions: {source_SVA_path}')
    print(f'File with VNTR motifs: {motifs_path}')
    print(f'File with reference genome: {reference_fasta_path}')
    print(f'File with SVA VNTR motifs: {SVA_VNTR_path}')
    
    # Get consensus sequences
    consensus_dict = consensus_seqs(consensus_path)
    # Open SVAs VNTR motifs file
    SVA_VNTR_motifs = read_file_and_store_lines(SVA_VNTR_path)
    # Get number or probabilities of events
    df_insertions1 = probabilities_total_number(probabilities_numbers_path, num_events)
    # Process insertion features and generate dict with random numbers based on distributions for each feature of every event
    dict_random = process_insertion_features_random_numbers(insertion_features_path, num_events, consensus_dict)
    # Add chromosome and start position
    df_insertions2 = add_beg_end_columns(df_insertions1, genome_wide_path)
    # Add the features of each event
    df_insertions3 = add_elements_columns(dict_random, df_insertions2)
    # Add SVA events additional details
    df_insertions4 = add_SVA_info(df_insertions3, consensus_dict)
    # Update df
    df_insertions5 = update_dataframe(df_insertions4, consensus_dict)
    # Add source gene information for transduction events
    df_insertions6 = add_source_gene_info(df_insertions5, source_L1_path, source_SVA_path)
    # Generate the insertion sequence and selected VNTR motifs
    df_insertions6['Sequence_Insertion'], df_insertions6['Selected_VNTR_Motifs'],df_insertions6['VNTR_Num_Motifs'] = zip(*df_insertions6.apply(lambda row: generate_insertion_seq(row, motifs_path, reference_fasta_path, consensus_dict, SVA_VNTR_motifs), axis=1))
    #df_insertions6['Sequence_Insertion'] = df_insertions6.apply(lambda row: generate_insertion_seq(row, motifs_path, reference_fasta_path, consensus_dict, SVA_VNTR_motifs), axis=1)
    # Update VNTR motifs
    df_insertions6 = process_vntr_motifs(df_insertions6)
    # Update the insertion sequence
    df_insertions6 = update_sequences(df_insertions6)
    # Save the output
    df_insertions6.to_csv('Insertions_table.tsv', sep='\t', index=False)
    
    # If the VCF argument is provided, create a VCF file
    if apply_VCF:
        df_insertions7 = pd.read_csv('Insertions_table.tsv', sep='\t')
        # Process the df to transform it to VCF format
        df_VCF_format = df_VCF(df_insertions7, reference_fasta_path)
        # Create the VCF
        create_vcf_file(df_VCF_format, reference_fasta_path, chromosome_length_path)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Generate insertion sequences returned in tsv file.')
    parser.add_argument('--consensus_path', type=str, required=True, help='Path to file with consensus sequences.')
    parser.add_argument('--probabilities_numbers_path', type=str, required=True, help='Path to the TSV file probabilities or defined number of events.')
    parser.add_argument('--insertion_features_path', type=str, required=True, help='Path to the TSV file containing insertion features of events.')
    parser.add_argument('--genome_wide_path', type=str, required=True, help='Path to the TSV file containing genome-wide distribution of events.')
    parser.add_argument('--source_L1_path', type=str, required=True, help='Path to the TSV file containing loci for LINE-1 transductions.')
    parser.add_argument('--source_SVA_path', type=str, required=True, help='Path to the TSV file containing loci for SVA transductions.')
    parser.add_argument('--motifs_path', type=str, required=True, help='Path to the TSV file containing loci for SVA transductions.')
    parser.add_argument('--SVA_VNTR_path', type=str, required=True, help='Path to the TSV file containing loci for SVA transductions.')
    parser.add_argument('--reference_fasta_path', type=str, required=True, help='Path to file with reference genome.')
    parser.add_argument('--chromosome_length_path', type=str, required=True, help='Path to the chromosome length file.')
    parser.add_argument('--num_events', type=int, default=100, required=False, help='Number of events to sample (optional, just in case of providing probabilities).')
    parser.add_argument('--VCF', action='store_true', required=False, help='If specified, creates a Variant Calling File (VCF)')
    
    args = parser.parse_args()

    main(args.consensus_path, args.probabilities_numbers_path, args.insertion_features_path, args.genome_wide_path, args.source_L1_path, args.source_SVA_path, args.motifs_path, args.SVA_VNTR_path, args.reference_fasta_path, args.chromosome_length_path, args.num_events, apply_VCF=args.VCF) 
import warnings

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
        dist.fit_transform(value, verbose = 0)
        random_numbers = dist.generate(num_samples)
        random_numbers = random_numbers.astype(int)
        random_numbers_dict[key] = random_numbers
    return random_numbers_dict

def remove_negative_values(values):
    ''' 
    Function to remove negative values from the list
    '''
    return [value for value in values if value > 0]

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
    Function to generate VNTR sequences by selecting a random row from the provided file
    and using its 'Complete_Sequence', 'Start', 'VNTR_Num_Motifs', and 'VNTR_Motifs'.
    '''

    # Step 1: Read the motifs file
    motifs_df = pd.read_csv(motifs_file, sep='\t')  # Reads the .tsv file containing motif information
    
    # Step 2: Randomly select a row from the motifs file
    random_row = motifs_df.sample(n=1).iloc[0]
    
    # Step 3: Get values from the randomly selected row
    complete_sequence = random_row['Complete_Sequence']
    start_position = random_row['Start']
    vntr_num_motifs = int(random_row['VNTR_Num_Motifs'])
    vntr_motifs = random_row['VNTR_Motifs'].split(',')
    # Step 4: Update the row with the selected values
    row['Sequence_Insertion'] = complete_sequence
    row['Length'] = len(complete_sequence)  # Set the Length as the length of the Complete_Sequence
    row['VNTR_Num_Motifs'] = vntr_num_motifs  # Update the VNTR_Num_Motifs
    row['VNTR_Motifs'] = vntr_motifs  # Update the VNTR_Motifs list
    row['Start'] = start_position  # Update the Start position
    
    # Step 5: Return the Complete_Sequence as the sequence and the VNTR_Motifs
    sequence = complete_sequence  # Now the sequence is just the Complete_Sequence
    return sequence, vntr_motifs, vntr_num_motifs

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
        rev_seq = reverse_complementary(rev_seq)  # Apply reverse complementary to REV sequence
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
        rev_seq = reverse_complementary(rev_seq)  # Apply reverse complementary to REV sequence
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
        rev_seq = reverse_complementary(rev_seq)  # Apply reverse complementary to REV sequence
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
        rev_seq = reverse_complementary(rev_seq)  # Apply reverse complementary to REV sequence
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
        rev_seq = reverse_complementary(rev_seq)  # Apply reverse complementary to REV sequence
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
        rev_seq = reverse_complementary(rev_seq)  # Apply reverse complementary to REV sequence
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
        sequence, selected_motifs, vntr_num_motifs = VNTR_insertions(row, motifs_file)
        return sequence, selected_motifs, vntr_num_motifs

    # DUPLICATIONS
    elif row['name'] == 'DUP':
        return DUP_insertions(row, reference_fasta), 0, 0
    
    # NUMT
    elif row['name'] == 'NUMT':
        return NUMT_insertions(row, dict_consensus), 0, 0

    # INVERSE DUPLICATIONS  
    elif row['name'] == 'INV_DUP':
        seq = DUP_insertions(row, reference_fasta)
        return reverse_complementary(seq), 0, 0

    # ORPHAN
    elif row['name'] == 'orphan':
        return orphan_insertions(row, reference_fasta)   , 0, 0

    # Alu__FOR+POLYA
    elif row['name'] == 'Alu__FOR+POLYA':
        return Alu__FOR_POLYA_insertions(row, dict_consensus, reference_fasta)   , 0, 0

    # Alu__TRUN+FOR+POLYA
    elif row['name'] == 'Alu__TRUN+FOR+POLYA':
        return Alu__FOR_POLYA_insertions(row, dict_consensus, reference_fasta)   , 0, 0

    # L1__FOR+POLYA
    elif row['name'] == 'L1__FOR+POLYA':
        return L1__FOR_POLYA_insertions(row, dict_consensus, reference_fasta)   , 0, 0
    
    # L1__TRUN+FOR+POLYA
    elif row['name'] == 'L1__TRUN+FOR+POLYA':
        return L1__FOR_POLYA_insertions(row, dict_consensus, reference_fasta) , 0  , 0
    
    # L1__TD+FOR+POLYA
    elif row['name'] == 'L1__TD+FOR+POLYA':
        return L1__TD_FOR_POLYA_insertions(row, dict_consensus, reference_fasta)   , 0, 0

    # L1__FOR+POLYA+TD+POLYA
    elif row['name'] == 'L1__FOR+POLYA+TD+POLYA':
        return L1__FOR_POLYA_TD_POLYA_insertions(row, dict_consensus, reference_fasta)   , 0, 0

    # L1__TRUN+FOR+POLYA+TD+POLYA
    elif row['name'] == 'L1__TRUN+FOR+POLYA+TD+POLYA':
        return L1__FOR_POLYA_TD_POLYA_insertions(row, dict_consensus, reference_fasta)   , 0, 0

    # L1__TRUN+REV+BLUNT+FOR+POLYA
    elif row['name'] == 'L1__TRUN+REV+BLUNT+FOR+POLYA':
        return L1__TRUN_REV_BLUNT_FOR_POLYA_insertions(row, dict_consensus, reference_fasta) , 0, 0

    # L1__TRUN+REV+BLUNT+FOR+POLYA+TD+POLYA
    elif row['name'] == 'L1__TRUN+REV+BLUNT+FOR+POLYA+TD+POLYA':
        return L1__TRUN_REV_BLUNT_FOR_POLYA_TD_POLYA_insertions(row, dict_consensus, reference_fasta) , 0, 0
    
    # L1__TRUN+REV+DUP+FOR+POLYA
    elif row['name'] == 'L1__TRUN+REV+DUP+FOR+POLYA':
        return L1__TRUN_REV_DUP_FOR_POLYA_insertions(row, dict_consensus, reference_fasta) , 0, 0

    # L1__TRUN+REV+DUP+FOR+POLYA+TD+POLYA
    elif row['name'] == 'L1__TRUN+REV+DUP+FOR+POLYA+TD+POLYA':
        return L1__TRUN_REV_DUP_FOR_POLYA_TD_POLYA_insertions(row, dict_consensus, reference_fasta) , 0, 0

    # L1__TRUN+REV+DEL+FOR+POLYA
    elif row['name'] == 'L1__TRUN+REV+DEL+FOR+POLYA':
        return L1__TRUN_REV_DEL_FOR_POLYA_insertions(row, dict_consensus, reference_fasta) , 0, 0
    
    # L1__REV+DEL+FOR+POLYA
    elif row['name'] == 'L1__REV+DEL+FOR+POLYA':
        return L1__TRUN_REV_DEL_FOR_POLYA_insertions(row, dict_consensus, reference_fasta) , 0, 0

    # L1__TRUN+REV+DEL+FOR+POLYA+TD+POLYA
    elif row['name'] == 'L1__TRUN+REV+DEL+FOR+POLYA+TD+POLYA':
        return L1__TRUN_REV_DEL_FOR_POLYA_TD_POLYA_insertions(row, dict_consensus, reference_fasta) , 0, 0

    # SVA__SINE-R+POLYA
    elif row['name'] == 'SVA__SINE-R+POLYA':
        return SVA__SINE_R_POLYA_insertions(row, dict_consensus, reference_fasta) , 0, 0

    # SVA__VNTR+SINE-R+POLYA
    elif row['name'] == 'SVA__VNTR+SINE-R+POLYA':
        return SVA__VNTR_SINE_R_POLYA_insertions(row, dict_consensus, reference_fasta, SVA_VNTR_list)   , 0 , 0

    # SVA__Alu-like+VNTR+SINE-R+POLYA
    elif row['name'] == 'SVA__Alu-like+VNTR+SINE-R+POLYA':
        return SVA__Alu_like_VNTR_SINE_R_POLYA_insertions(row, dict_consensus, reference_fasta, SVA_VNTR_list)    , 0, 0

    # SVA__MAST2+VNTR+SINE-R+POLYA
    elif row['name'] == 'SVA__MAST2+VNTR+SINE-R+POLYA':
        return SVA__MAST2_VNTR_SINE_R_POLYA_insertions(row, dict_consensus, reference_fasta, SVA_VNTR_list)  , 0, 0

    # SVA__TD+MAST2+VNTR+SINE-R+POLYA
    elif row['name'] == 'SVA__TD+MAST2+VNTR+SINE-R+POLYA':
        return SVA__TD_MAST2_VNTR_SINE_R_POLYA_insertions(row, dict_consensus, reference_fasta, SVA_VNTR_list)  , 0  , 0

    # SVA__SINE-R+POLYA+TD+POLYA
    elif row['name'] == 'SVA__SINE-R+POLYA+TD+POLYA':
        return SVA__SINE_R_POLYA_TD_POLYA_insertions(row, dict_consensus, reference_fasta)   , 0, 0

    # SVA__VNTR+SINE-R+POLYA+TD+POLYA
    elif row['name'] == 'SVA__VNTR+SINE-R+POLYA+TD+POLYA':
        return SVA__VNTR_SINE_R_POLYA_TD_POLYA_insertions(row, dict_consensus, reference_fasta, SVA_VNTR_list), 0, 0

    # SVA__Alu-like+VNTR+SINE-R+POLYA+TD+POLYA
    elif row['name'] == 'SVA__Alu-like+VNTR+SINE-R+POLYA+TD+POLYA':
        return SVA__Alu_like_VNTR_SINE_R_POLYA_TD_POLYA_insertions(row, dict_consensus, reference_fasta, SVA_VNTR_list), 0, 0

    # SVA__MAST2+VNTR+SINE-R+POLYA+TD+POLYA
    elif row['name'] == 'SVA__MAST2+VNTR+SINE-R+POLYA+TD+POLYA':
        return SVA__MAST2_VNTR_SINE_R_POLYA_TD_POLYA_insertions(row, dict_consensus, reference_fasta, SVA_VNTR_list)  , 0, 0

    # SVA__Hexamer+Alu-like+VNTR+SINE-R+POLYA
    elif row['name'] == 'SVA__Hexamer+Alu-like+VNTR+SINE-R+POLYA':
        return SVA__Hexamer_Alu_like_VNTR_SINE_R_POLYA_insertions(row, dict_consensus, reference_fasta, SVA_VNTR_list)    , 0, 0

    # SVA__TD+Hexamer+Alu-like+VNTR+SINE-R+POLYA
    elif row['name'] == 'SVA__TD+Hexamer+Alu-like+VNTR+SINE-R+POLYA':
        return SVA__TD_Hexamer_Alu_like_VNTR_SINE_R_POLYA_insertions(row, dict_consensus, reference_fasta, SVA_VNTR_list)    , 0, 0
    
    # SVA__Hexamer+Alu-like+VNTR+SINE-R+POLYA+TD+POLYA
    elif row['name'] == 'SVA__Hexamer+Alu-like+VNTR+SINE-R+POLYA+TD+POLYA':
        return SVA__Hexamer_Alu_like_VNTR_SINE_R_POLYA_TD_POLYA_insertions(row, dict_consensus, reference_fasta, SVA_VNTR_list)   , 0 , 0

    else:
        return '', 0, 0 

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

# Function to process each row's 'name' value
def process_name(name):
    if '__' in name:  # Case when there is '__' in the name
        fam_n, conformation = name.split('__', 1)
        itype_n = 'partnered' if 'TD' in conformation else 'solo'
    else:  # Case when there is no '__' in the name
        fam_n, conformation = np.nan, np.nan
        itype_n = name
    
    return pd.Series([itype_n, conformation, fam_n])

# Function to generate the HEXAMER_SEQ column
def generate_hexamer_seq(row):
    # Check if 'FAM_N' contains 'SVA' and 'CONFORMATION' contains 'Hexamer'
    if pd.notna(row['FAM_N']) and 'SVA' in row['FAM_N'] and pd.notna(row['CONFORMATION']) and 'Hexamer' in row['CONFORMATION']:
        hexamer = 'CCCTCT'
        try:
            # Assuming HEXAMER_LEN is a column with the desired hexamer length
            sva_hexamer_length = int(row['HEXAMER_LEN'])  # Desired length of hexamer
            # Generate the hexamer sequence
            hexamer_seq = hexamer * (sva_hexamer_length // 6) + hexamer[:sva_hexamer_length % 6]
            return hexamer_seq
        except (ValueError, TypeError):
            return np.nan  # Return NaN if there's an issue with the length
    return np.nan  # Return NaN if conditions are not met

# Function to generate TSD sequence from the reference fasta file
def generate_tsd_seq(row, reference_fasta):
    # Check if 'beg' and 'TSD_LEN' are valid
    if pd.notna(row['beg']) and pd.notna(row['TSD_LEN']):
        try:
            beg = int(row['beg'])
            tsd_len = int(row['TSD_LEN'])
            start_tsd = beg - tsd_len

            # Fetch the sequence from the reference FASTA file using pysam
            with pysam.FastaFile(reference_fasta) as fasta_file:
                TSD = fasta_file.fetch(row['#ref'], start_tsd, beg)  # Extract the sequence
                return TSD.upper()
        except (ValueError, TypeError, KeyError):
            return np.nan  # Return NaN if there's an issue with the length or values
    return np.nan  # Return NaN if either 'beg' or 'TSD_LEN' is missing or invalid

# Function to generate sequence for 3PRIME_TD
def generate_3prime_td_seq(row, reference_fasta):
    # Fetch the sequence for 3PRIME_TD_SEQ if there's a value in 3PRIME_NB_TD
    if pd.notna(row['3PRIME_NB_TD']):
        try:
            with pysam.FastaFile(reference_fasta) as fasta_file:
                transduction = fasta_file.fetch(row['SRC_ref'], row['TD_beg'], row['TD_end'])
                return transduction
        except (KeyError, ValueError):
            return np.nan
    return np.nan

# Function to generate sequence for 5PRIME_TD
def generate_5prime_td_seq(row, reference_fasta):
    # Fetch the sequence for 5PRIME_TD_SEQ if there's a value in 5PRIME_NB_TD
    if pd.notna(row['5PRIME_NB_TD']):
        try:
            with pysam.FastaFile(reference_fasta) as fasta_file:
                transduction = fasta_file.fetch(row['SRC_ref'], row['TD_beg'], row['TD_end'])
                return transduction
        except (KeyError, ValueError):
            return np.nan
    return np.nan

# Function to generate sequence for orphan TD
def generate_orphan_seq(row, reference_fasta):
    # Fetch the sequence for orphan TD if there's a value in ORPHAN_TD_COORD
    if pd.notna(row['ORPHAN_TD_COORD']):
        try:
            with pysam.FastaFile(reference_fasta) as fasta_file:
                transduction = fasta_file.fetch(row['SRC_ref'], row['TD_beg'], row['TD_end'])
                return transduction
        except (KeyError, ValueError):
            return np.nan
    return np.nan

# Function to generate POLYA_LEN column
def generate_polya_len(row):
    # Remove the decimal by converting the value to an integer before adding to the string
    poly_a_len_1 = str(int(row['PolyA_Length_1'])) if pd.notna(row['PolyA_Length_1']) else ''
    poly_a_len_2 = str(int(row['PolyA_Length_2'])) if pd.notna(row['PolyA_Length_2']) else ''
    
    if poly_a_len_1 and poly_a_len_2:
        return poly_a_len_1 + ',' + poly_a_len_2
    elif poly_a_len_1:
        return poly_a_len_1
    else:
        return np.nan  # Return NaN if both are missing

# Function to generate POLYA_SEQ column
def generate_polya_seq(row):
    strand = row['STRAND'] if pd.notna(row['STRAND']) else '+'
    
    # Determine the PolyA sequences
    seq_1 = 'A' * int(row['PolyA_Length_1']) if pd.notna(row['PolyA_Length_1']) else ''
    seq_2 = 'A' * int(row['PolyA_Length_2']) if pd.notna(row['PolyA_Length_2']) else ''
    
    if strand == '-':
        seq_1 = seq_1.replace('A', 'T') if seq_1 else ''
        seq_2 = seq_2.replace('A', 'T') if seq_2 else ''
    
    # Join the sequences with a comma if both exist
    if seq_1 and seq_2:
        return seq_1 + ',' + seq_2
    elif seq_1:
        return seq_1
    elif seq_2:
        return seq_2
    return np.nan  # Return NaN if no sequences are available

def process_vntr_motifs(df):
    # Iterate over each element in the 'Selected_VNTR_Motifs' column
    for i in range(len(df)):
        #motif = df.at[i, 'Selected_VNTR_Motifs']
        motif = df.iloc[i]['Selected_VNTR_Motifs']
        if motif == '0':
            df.at[i, 'Selected_VNTR_Motifs'] = np.nan  # Replace '0' with NaN
        elif isinstance(motif, str):
            # Remove '[' and ']' from the string
            df.at[i, 'Selected_VNTR_Motifs'] = motif.replace("[", "").replace("'[", "").replace("]", "").replace("]'", "").strip()
    
    return df

# Main function
def df_VCF(df, reference_fasta):
    # Create a copy of the DataFrame
    df_copy = df.copy()
    df_copy = process_vntr_motifs(df_copy)
    df_copy['Selected_VNTR_Motifs'] = df_copy['Selected_VNTR_Motifs'].str.replace("'", "")
    # Rename columns based on your requirements
    df_copy = df_copy.rename(columns={
        'Length': 'INS_LEN',
        'VNTR_Num_Motifs': 'NB_MOTIFS',
        'SVA_Hexamer': 'HEXAMER_LEN',
        'TSD_Length': 'TSD_LEN',
        'TD_5': '5PRIME_TD_LEN',
        'TD_3': '3PRIME_TD_LEN',
        'Strand': 'STRAND',
        'TD_orphan_Length': 'ORPHAN_TD_LEN',
        'Selected_VNTR_Motifs': 'MOTIFS'
    })

    # Convert columns to integers to remove decimals (i.e., .0)
    df_copy['NB_MOTIFS'] = df_copy['NB_MOTIFS'].apply(lambda x: int(x) if pd.notna(x) else np.nan).astype('Int64')  # Use 'Int64' to keep NaNs
    df_copy['HEXAMER_LEN'] = df_copy['HEXAMER_LEN'].apply(lambda x: int(x) if pd.notna(x) else np.nan).astype('Int64')  # Use 'Int64' to keep NaNs
    df_copy['TSD_LEN'] = df_copy['TSD_LEN'].apply(lambda x: int(x) if pd.notna(x) else np.nan).astype('Int64')  # Use 'Int64' to keep NaNs
    df_copy['5PRIME_TD_LEN'] = df_copy['5PRIME_TD_LEN'].apply(lambda x: int(x) if pd.notna(x) else np.nan).astype('Int64')  # Use 'Int64' to keep NaNs
    df_copy['3PRIME_TD_LEN'] = df_copy['3PRIME_TD_LEN'].apply(lambda x: int(x) if pd.notna(x) else np.nan).astype('Int64')  # Use 'Int64' to keep NaNs
    
    # Add the 'ID' column with the format SV1, SV2, SV3...
    df_copy['ID'] = ['INS_' + str(i + 1) for i in range(len(df_copy))]

    # Apply the processing function to the 'name' column and split it into three new columns
    df_copy[['ITYPE_N', 'CONFORMATION', 'FAM_N']] = df_copy['name'].apply(process_name)
    
    # For rows where there is no 'CONFORMATION' (e.g., VNTR), set 'CONFORMATION' to NaN
    df_copy['CONFORMATION'] = df_copy['CONFORMATION'].replace('', np.nan)

    # Remove the 'name' column at the end
    df_copy = df_copy.drop(columns=['name'])
    
    # Apply the hexamer sequence generation logic
    df_copy['HEXAMER_SEQ'] = df_copy.apply(generate_hexamer_seq, axis=1)
    
    # Apply the TSD sequence generation logic
    df_copy['TSD_SEQ'] = df_copy.apply(lambda row: generate_tsd_seq(row, reference_fasta), axis=1)
    
    # Create the 3PRIME_NB_TD and 5PRIME_NB_TD columns based on the conditions
    df_copy['3PRIME_NB_TD'] = df_copy['3PRIME_TD_LEN'].apply(lambda x: 1 if pd.notna(x) else np.nan)
    df_copy['5PRIME_NB_TD'] = df_copy['5PRIME_TD_LEN'].apply(lambda x: 1 if pd.notna(x) else np.nan)
    
    # Create the 3PRIME_TD_COORD, 5PRIME_TD_COORD and ORPHAN_TD_COORD columns based on the conditions
    df_copy['3PRIME_TD_COORD'] = df_copy['SRC_identifier'].where(df_copy['3PRIME_NB_TD'].notna(), np.nan)
    df_copy['5PRIME_TD_COORD'] = df_copy['SRC_identifier'].where(df_copy['5PRIME_NB_TD'].notna(), np.nan)
    df_copy['ORPHAN_TD_COORD'] = df_copy['SRC_identifier'].where(df_copy['ITYPE_N'] == 'orphan')

    # Remove the SRC_identifier column
    df_copy = df_copy.drop(columns=['SRC_identifier'])

    # Apply the transduction sequence generation logic for 3PRIME_TD_SEQ and 5PRIME_TD_SEQ
    df_copy['3PRIME_TD_SEQ'] = df_copy.apply(lambda row: generate_3prime_td_seq(row, reference_fasta), axis=1)
    df_copy['5PRIME_TD_SEQ'] = df_copy.apply(lambda row: generate_5prime_td_seq(row, reference_fasta), axis=1)
    df_copy['ORPHAN_TD_SEQ'] = df_copy.apply(lambda row: generate_orphan_seq(row, reference_fasta), axis=1)  

    # Convert the 3PRIME_TD_SEQ and 5PRIME_TD_SEQ columns to uppercase
    df_copy['3PRIME_TD_SEQ'] = df_copy['3PRIME_TD_SEQ'].apply(lambda x: x.upper() if pd.notna(x) else np.nan)
    df_copy['5PRIME_TD_SEQ'] = df_copy['5PRIME_TD_SEQ'].apply(lambda x: x.upper() if pd.notna(x) else np.nan)
    df_copy['ORPHAN_TD_SEQ'] = df_copy['ORPHAN_TD_SEQ'].apply(lambda x: x.upper() if pd.notna(x) else np.nan)

    # Apply the POLYA_LEN generation logic
    df_copy['POLYA_LEN'] = df_copy.apply(generate_polya_len, axis=1)

    # Apply the POLYA_SEQ generation logic
    df_copy['POLYA_SEQ'] = df_copy.apply(generate_polya_seq, axis=1)

    return df_copy

def create_vcf_file(df, reference_fasta, chromosome_length):
    # Create the chr_length dictionary using formats
    chr_length = formats.chrom_lengths_index(chromosome_length)

    # Get date of creation
    current_date = datetime.datetime.now().strftime("%Y%m%d")

    # Open a VCF file to write to
    with open('VCF_Insertions_SVModeller.vcf', 'w') as vcf_file:
        # Write VCF header
        vcf_file.write("##fileformat=VCFv4.2\n")
        vcf_file.write(f"##fileDate={current_date}\n")
        vcf_file.write("##source=SVModeller\n")
        vcf_file.write("##reference=None\n")

        # Adding contig length
        for contig in df['#ref'].unique():
            # Get the length from chr_length dictionary
            if contig in chr_length:
                contig_length = chr_length[contig]
                vcf_file.write(f"##contig=<ID={contig},assembly=None,length={contig_length},species=None>\n")
            else:
                # Handle the case where the contig is not found in chr_length
                vcf_file.write(f"##contig=<ID={contig},assembly=None,length=None,species=None>\n")

        vcf_file.write("##INFO=<ID=INS_LEN,Type=float,Description=Total length of the insertion>\n")
        vcf_file.write("##INFO=<ID=TSD_LEN,Type=float,Description=Length of the TSD>\n")
        vcf_file.write("##INFO=<ID=5PRIME_TD_LEN,Type=float,Description=Length of the 5' TD>\n")
        vcf_file.write("##INFO=<ID=3PRIME_TD_LEN,Type=float,Description=Length of the 3' TD>\n")
        vcf_file.write("##INFO=<ID=NB_MOTIFS,Type=float,Description=Number of VNTR motifs>\n")
        vcf_file.write("##INFO=<ID=MOTIFS,Type=str,Description=VNTR selected motifs>\n")
        vcf_file.write("##INFO=<ID=HEXAMER_LEN,Type=float,Description=Length of the SVA hexamer>\n")
        vcf_file.write("##INFO=<ID=SVA_VNTR_Length,Type=float,Description=Length of the SVA VNTR>\n")
        vcf_file.write("##INFO=<ID=ITYPE_N,Type=float,Description=Type of insertion>\n")
        vcf_file.write("##INFO=<ID=CONFORMATION,Type=float,Description=Conformation of the insertion>\n")
        vcf_file.write("##INFO=<ID=FAM_N,Type=float,Description=Family of the insertion>\n")
        vcf_file.write("##INFO=<ID=HEXAMER_SEQ,Type=float,Description=Sequence of SVA hexamer>\n")
        vcf_file.write("##INFO=<ID=TSD_SEQ,Type=float,Description=Length of the poly A tail (or first poly A in case of more than 1)>\n")
        vcf_file.write("##INFO=<ID=3PRIME_NB_TD,Type=float,Description=Number of 3' TD>\n")
        vcf_file.write("##INFO=<ID=5PRIME_NB_TD,Type=float,Description=Number of 5' TD>\n")
        vcf_file.write("##INFO=<ID=3PRIME_TD_COORD,Type=float,Description=Coordinates of 3' TD>\n")
        vcf_file.write("##INFO=<ID=5PRIME_TD_COORD,Type=float,Description=Coordinates of 5' TD>\n")
        vcf_file.write("##INFO=<ID=3PRIME_TD_SEQ,Type=float,Description=Sequence of 3' TD>\n")
        vcf_file.write("##INFO=<ID=5PRIME_TD_SEQ,Type=float,Description=Sequence of 5' TD>\n")
        vcf_file.write("##INFO=<ID=POLYA_LEN,Type=float,Description=Length of the poly A tail>\n")
        vcf_file.write("##INFO=<ID=POLYA_SEQ,Type=float,Description=Sequence of the poly A tail>\n")
        vcf_file.write("##INFO=<ID=FOR,Type=float,Description=Length of conformation forward part of the event>\n")
        vcf_file.write("##INFO=<ID=TRUN,Type=float,Description=Length of the conformation truncated part of the event>\n")
        vcf_file.write("##INFO=<ID=REV,Type=float,Description=Length of the conformation reverse part of the event>\n")
        vcf_file.write("##INFO=<ID=DEL,Type=float,Description=Length of the conformation deleted part of the event>\n")
        vcf_file.write("##INFO=<ID=DUP,Type=float,Description=Length of the conformation duplicated part of the event>\n")
        vcf_file.write("##INFO=<ID=STRAND,Type=str,Description=Strand (+ or -) of the event>\n")
        vcf_file.write("##INFO=<ID=ORPHAN_TD_LEN,Type=float,Description=Length of the orphan transduction>\n")
        vcf_file.write("##INFO=<ID=ORPHAN_TD_COORD,Type=str,Description=Coordinates of orphan transduction>\n")
        vcf_file.write("##INFO=<ID=ORPHAN_TD_SEQ,Type=str,Description=Sequence of orphan transduction>\n")

        vcf_file.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n")

        # Loop through each row in the DataFrame
        for _, row in df.iterrows():
            chrom = row['#ref']
            pos = row['beg']
            event_id = row['ID']
            ref = '.' 
            alt = row['Sequence_Insertion']
            qual = '.'  
            filter = '.'

            # Convert 'beg' to an integer (position at which to fetch the reference sequence)
            beg = int(row['beg'])

            # Fetch the sequence using pysam from the reference genome at the position 'beg'
            with pysam.FastaFile(reference_fasta) as fasta_file:
                ref = fasta_file.fetch(chrom, beg - 1, beg)  # pysam is 0-based, so we subtract 1 for 0-based indexing
            ref = ref.upper()

            # Create the INFO field dynamically, excluding NaN values
            info_fields = []
            for col in ['ITYPE_N', 'INS_LEN', 'STRAND', 'FAM_N', 'CONFORMATION', 'FOR', 'TRUN', 'REV', 'DEL', 'DUP', 'TSD_LEN', 'TSD_SEQ', '3PRIME_NB_TD', '5PRIME_NB_TD', '5PRIME_TD_LEN', '3PRIME_TD_LEN', '5PRIME_TD_COORD', '3PRIME_TD_COORD','3PRIME_TD_SEQ', '5PRIME_TD_SEQ', 'NB_MOTIFS', 
                        'MOTIFS','HEXAMER_LEN', 'SVA_VNTR_Length', 'HEXAMER_SEQ', 'ORPHAN_TD_LEN','ORPHAN_TD_COORD', 'ORPHAN_TD_SEQ','POLYA_LEN', 'POLYA_SEQ']:
                value = row[col]
                if pd.notna(value):  # Check if the value is not NaN
                    info_fields.append(f"{col}={value}")

            # Join the info fields into a single string separated by semicolons
            info = ";".join(info_fields)

            # Write the VCF entry for each row
            vcf_file.write(f"{chrom}\t{pos}\t{event_id}\t{ref}\t{alt}\t{qual}\t{filter}\t{info}\n")

    print("VCF file created successfully.")

# Remove FutureWarnings
warnings.simplefilter(action='ignore', category=FutureWarning)

def main(consensus_path, probabilities_numbers_path, insertion_features_path, genome_wide_path, source_L1_path, source_SVA_path, motifs_path, SVA_VNTR_path, reference_fasta_path, chromosome_length_path, num_events, apply_VCF):
    print(f'File with consensus sequences: {consensus_path}')
    print(f'File with probabilities or number of events: {probabilities_numbers_path}')
    print(f'File with insertions features: {insertion_features_path}')
    print(f'File with genome-wide distribution: {genome_wide_path}')
    print(f'File with source genes for L1 transductions: {source_L1_path}')
    print(f'File with source genes for SVA transductions: {source_SVA_path}')
    print(f'File with VNTR motifs: {motifs_path}')
    print(f'File with reference genome: {reference_fasta_path}')
    print(f'File with SVA VNTR motifs: {SVA_VNTR_path}')
    
    # Get consensus sequences
    consensus_dict = consensus_seqs(consensus_path)
    # Open SVAs VNTR motifs file
    SVA_VNTR_motifs = read_file_and_store_lines(SVA_VNTR_path)
    # Get number or probabilities of events
    df_insertions1 = probabilities_total_number(probabilities_numbers_path, num_events)
    # Process insertion features and generate dict with random numbers based on distributions for each feature of every event
    dict_random = process_insertion_features_random_numbers(insertion_features_path, num_events, consensus_dict)
    # Add chromosome and start position
    df_insertions2 = add_beg_end_columns(df_insertions1, genome_wide_path)
    # Add the features of each event
    df_insertions3 = add_elements_columns(dict_random, df_insertions2)
    # Add SVA events additional details
    df_insertions4 = add_SVA_info(df_insertions3, consensus_dict)
    # Update df
    df_insertions5 = update_dataframe(df_insertions4, consensus_dict)
    # Add source gene information for transduction events
    df_insertions6 = add_source_gene_info(df_insertions5, source_L1_path, source_SVA_path)
    # Generate the insertion sequence and selected VNTR motifs
    df_insertions6['Sequence_Insertion'], df_insertions6['Selected_VNTR_Motifs'],df_insertions6['VNTR_Num_Motifs'] = zip(*df_insertions6.apply(lambda row: generate_insertion_seq(row, motifs_path, reference_fasta_path, consensus_dict, SVA_VNTR_motifs), axis=1))
    #df_insertions6['Sequence_Insertion'] = df_insertions6.apply(lambda row: generate_insertion_seq(row, motifs_path, reference_fasta_path, consensus_dict, SVA_VNTR_motifs), axis=1)
    # Update VNTR motifs
    df_insertions6 = process_vntr_motifs(df_insertions6)
    # Update the insertion sequence
    df_insertions6 = update_sequences(df_insertions6)
    # Save the output
    df_insertions6.to_csv('Insertions_table.tsv', sep='\t', index=False)
    
    # If the VCF argument is provided, create a VCF file
    if apply_VCF:
        df_insertions7 = pd.read_csv('Insertions_table.tsv', sep='\t')
        # Process the df to transform it to VCF format
        df_VCF_format = df_VCF(df_insertions7, reference_fasta_path)
        # Create the VCF
        create_vcf_file(df_VCF_format, reference_fasta_path, chromosome_length_path)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Generate insertion sequences returned in tsv file.')
    parser.add_argument('--consensus_path', type=str, required=True, help='Path to file with consensus sequences.')
    parser.add_argument('--probabilities_numbers_path', type=str, required=True, help='Path to the TSV file probabilities or defined number of events.')
    parser.add_argument('--insertion_features_path', type=str, required=True, help='Path to the TSV file containing insertion features of events.')
    parser.add_argument('--genome_wide_path', type=str, required=True, help='Path to the TSV file containing genome-wide distribution of events.')
    parser.add_argument('--source_L1_path', type=str, required=True, help='Path to the TSV file containing loci for LINE-1 transductions.')
    parser.add_argument('--source_SVA_path', type=str, required=True, help='Path to the TSV file containing loci for SVA transductions.')
    parser.add_argument('--motifs_path', type=str, required=True, help='Path to the TSV file containing loci for SVA transductions.')
    parser.add_argument('--SVA_VNTR_path', type=str, required=True, help='Path to the TSV file containing loci for SVA transductions.')
    parser.add_argument('--reference_fasta_path', type=str, required=True, help='Path to file with reference genome.')
    parser.add_argument('--chromosome_length_path', type=str, required=True, help='Path to the chromosome length file.')
    parser.add_argument('--num_events', type=int, default=100, required=False, help='Number of events to sample (optional, just in case of providing probabilities).')
    parser.add_argument('--VCF', action='store_true', required=False, help='If specified, creates a Variant Calling File (VCF)')
    
    args = parser.parse_args()

    main(args.consensus_path, args.probabilities_numbers_path, args.insertion_features_path, args.genome_wide_path, args.source_L1_path, args.source_SVA_path, args.motifs_path, args.SVA_VNTR_path, args.reference_fasta_path, args.chromosome_length_path, args.num_events, apply_VCF=args.VCF)
