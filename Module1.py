# SVModeler

# Module 1.1

import argparse
import formats
import gRanges
import structures
import numpy as np
import statistics
import pandas as pd
import pysam
import math
import distfit
from distfit import distfit
import random
import matplotlib.pyplot as plt

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
    vcf = formats.VCF()  # Instantiate the VCF class
    vcf.read(file_path)  # Read the VCF file

    data = []  # List to store the extracted data
    
    for variant in vcf.variants:
        length = len(variant.alt)
        end_position = variant.pos + length
        
        # Not in all lines common variables
        # Get the x value, if not present, set it to NA

        poly_a_len = variant.info.get('POLYA_LEN', 'NA')  # Get the POLYA_LEN value, if not present, set it to 'NA'
        family = variant.info.get('FAM_N', 'NA')
        strand = variant.info.get('STRAND', 'NA')
        rt_len = variant.info.get('RT_LEN', 'NA')
        tsd_len = variant.info.get('TSD_LEN','NA')
        TD_3_Num = variant.info.get('3PRIME_NB_TD','NA')
        TD_3 = variant.info.get('3PRIME_TD_LEN','NA')
        TD_5_Num = variant.info.get('5PRIME_NB_TD','NA')
        TD_5 = variant.info.get('5PRIME_TD_LEN','NA')
        hexamer = variant.info.get('HEXAMER_LEN','NA')
        TD_orphan = variant.info.get('ORPHAN_TD_LEN','NA')
        VNTR_motif_num = variant.info.get('NB_MOTIFS','NA')
        VNTR_motif = variant.info.get('MOTIFS','NA')

        # Add to the data list the information
        data.append({
            'Type_SV': variant.info['ITYPE_N'], 
            'Family':family,
            'Start_position': variant.pos,
            'End_position': end_position,
            'Length': length,
            'Chromosome': variant.chrom,
            'PolyA_Length': poly_a_len,
            'Strand': strand,
            'RT_Length': rt_len,
            'TSD_Length': tsd_len,
            'TD_5_Num': TD_5_Num,
            'TD_5': TD_5,
            'TD_3_Num': TD_3_Num,
            'TD_3': TD_3,
            'SVA_Hexamer': hexamer,
            'TD_orphan_Length': TD_orphan,
            'VNTR_Num_Motifs': VNTR_motif_num,
            'VNTR_Motifs': VNTR_motif
        })
    
    # Create a data frame from the extracted data
    df = pd.DataFrame(data, columns=['Chromosome', 'Start_position', 'End_position', 'Type_SV', 'Family', 'Length', 'PolyA_Length', 'Strand', 'RT_Length', 'TSD_Length', 
                                     'TD_5_Num', 'TD_5', 'TD_3_Num', 'TD_3', 'SVA_Hexamer', 'TD_orphan_Length', 'VNTR_Num_Motifs', 'VNTR_Motifs'])
    
    return df


def TD_filter(df):
    ''' 
    Input & output: df
    1- checks if TD_5_Num is not NA and has a value higher than 1
    2- if so, deletes the value of TD_5_Num and TD_5
    3- does the same for TD_3_Num
    '''
    for index, row in df.iterrows():
        if row['TD_5_Num'] != 'NA' and int(row['TD_5_Num']) > 1:
            df.at[index, 'TD_5'] = 'NA'
            df.at[index, 'TD_5_Num'] = 'NA'
        
        if row['TD_3_Num'] != 'NA' and int(row['TD_3_Num']) > 1:
            df.at[index, 'TD_3'] = 'NA'
            df.at[index, 'TD_3_Num'] = 'NA'
    return df 




def dict_mutations_list(df):
    ''' 
    Input: df in BED format
    Output: dictionary

    Function that takes a df in BED format and checks for each row if there are values of the specific columns and adds those
    values to a new dictionary that has as keys the type of SV + PolyA,Length,...
    '''
    # Create an empty dictionary
    info_dict = {}

    # Iterate over all rows of the dataframe
    for index, row in df.iterrows():
        if row['PolyA_Length_1'] != 'NA':
            key = row['name'] + '_' + 'PolyA_Length_1'
            if key in info_dict:
                info_dict[key].append(row['PolyA_Length_1'])
            else:
                info_dict[key] = [row['PolyA_Length_1']]

        if row['PolyA_Length_2'] != 'NA':
            key = row['name'] + '_' + 'PolyA_Length_2'
            if key in info_dict:
                info_dict[key].append(row['PolyA_Length_2'])
            else:
                info_dict[key] = [row['PolyA_Length_2']]

        if row['RT'] != ('NA' or 'nan'):
            key = row['name'] + '_' + 'RT'
            if key in info_dict:
                info_dict[key].append(row['RT'])
            else:
                info_dict[key] = [row['RT']]
        
        if row['TSD_Length'] != 'NA':
            key = row['name'] + '_' + 'TSD_Length'
            if key in info_dict:
                info_dict[key].append(row['TSD_Length'])
            else:
                info_dict[key] = [row['TSD_Length']]
        
        if row['Length'] != 'NA':
            key = row['name'] + '_' + 'Length'
            if key in info_dict:
                info_dict[key].append(row['Length'])
            else:
                info_dict[key] = [row['Length']]

        if row['TD_5'] != 'NA':
            key = row['name'] + '_' + 'TD_5'
            if key in info_dict:
                info_dict[key].append(row['TD_5'])
            else:
                info_dict[key] = [row['TD_5']]
            
        if row['TD_3'] != 'NA':
            key = row['name'] + '_' + 'TD_3'
            if key in info_dict:
                info_dict[key].append(row['TD_3'])
            else:
                info_dict[key] = [row['TD_3']]

        if row['SVA_Hexamer'] != 'NA':
            key = row['name'] + '_' + 'SVA_Hexamer'
            if key in info_dict:
                info_dict[key].append(row['SVA_Hexamer'])
            else:
                info_dict[key] = [row['SVA_Hexamer']]
        
        if row['TD_orphan_Length'] != 'NA':
            key = row['name'] + '_' + 'TD_orphan_Length'
            if key in info_dict:
                info_dict[key].append(row['TD_orphan_Length'])
            else:
                info_dict[key] = [row['TD_orphan_Length']]

        if row['VNTR_Num_Motifs'] != 'NA':
            key = row['name'] + '_' + 'VNTR_Num_Motifs'
            if key in info_dict:
                info_dict[key].append(row['VNTR_Num_Motifs'])
            else:
                info_dict[key] = [row['VNTR_Num_Motifs']]

        if row['VNTR_Motifs'] != 'NA':
            key = row['name'] + '_' + 'VNTR_Motifs'
            if key in info_dict:
                info_dict[key].append(row['VNTR_Motifs'])
            else:
                info_dict[key] = [row['VNTR_Motifs']]
                    
    return(info_dict)


def filter_sd(dict,key_list):
    ''' 
    Function to remove those values that exceed twice the standard deviation
    '''
    for key in key_list:
        values = dict[key]
        std_dev = statistics.stdev(values)
        values = [value for value in values if value <= 2 * std_dev]
        dict[key] = values


def remove_nan(lst):
    ''' 
    Function to remove nan values from a list
    '''
    return [x for x in lst if not np.isnan(x)]


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
        'partnered_Alu': 'partnered_Alu',
        'partnered_SVA': 'partnered_SVA',
        'partnered_L1': 'partnered_L1',
        'orphan': 'orphan',
        'PSD': 'PSD',
        'chimera': 'chimera',
        'partnered': 'partnered',
        'VNTR': 'VNTR',
        'DUP': 'DUP',
        'INV_DUP': 'INV_DUP',
        'COMPLEX_DUP': 'COMPLEX_DUP',
        'NUMT': 'NUMT'
    }

    # Initialize the output DataFrame
    df = pd.DataFrame(columns=['window', 'beg', 'end'] + list(subtype_cols.values()))

    # Iterate over each window
    for window in bins:
        chrom, start, end = window
        window_df = result_BED_table[(result_BED_table['#ref'] == chrom) & (result_BED_table['beg'] >= start) & (result_BED_table['end'] <= end)]

        # Count the occurrences of each SubType in the window
        counts = window_df['name'].value_counts().to_dict()

        # Create a row for the output DataFrame
        row = [chrom, start, end]
        for subtype in subtype_cols:
            row.append(counts.get(subtype, 0))

        # Append the row to the output DataFrame
        df.loc[len(df)] = row

    return df


def calculate_normalized_values(res_table):
    ''' 
    Function to calculate the normalized values of each variant of each chromosome
    '''
    # Get unique values in the 'window' column
    unique_chr_values = res_table['window'].unique()
    
    # Initialize a new DataFrame to store the results
    new_table = pd.DataFrame(columns=res_table.columns)
    
    # Iterate over each unique 'window' value
    for chr_value in unique_chr_values:
        # Filter the rows for the current 'window' value
        subset = res_table[res_table['window'] == chr_value]
        
        # Copy the values of 'beg' and 'end' columns
        normalized_values = subset[['beg', 'end']]
        
        # Calculate the sum of values for each column
        total_values = subset[['solo_SVA', 'solo_Alu', 'partnered_Alu', 'partnered_SVA',
       'partnered_L1', 'orphan', 'solo_L1', 'PSD', 'chimera', 'partnered',
       'VNTR', 'DUP', 'INV_DUP', 'COMPLEX_DUP', 'NUMT']].sum()
        
        # Divide each cell value by the total sum
        normalized_values[['solo_SVA', 'solo_Alu', 'partnered_Alu', 'partnered_SVA',
       'partnered_L1', 'orphan', 'solo_L1', 'PSD', 'chimera', 'partnered',
       'VNTR', 'DUP', 'INV_DUP', 'COMPLEX_DUP', 'NUMT']] = subset[['solo_SVA', 'solo_Alu', 'partnered_Alu', 'partnered_SVA',
       'partnered_L1', 'orphan', 'solo_L1', 'PSD', 'chimera', 'partnered',
       'VNTR', 'DUP', 'INV_DUP', 'COMPLEX_DUP', 'NUMT']].div(total_values).values
        
        # Add the 'window' values back to the DataFrame
        normalized_values['window'] = chr_value
        
        # Append the normalized values to the new DataFrame
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


def add_beg_end_columns(df1, df2):
    ''' 
    Function to add start and end columns to the df based on probabilities
    '''

    # Create new columns in the first DataFrame
    df1['#ref'] = ''
    df1['beg'] = ''
    df1['end'] = ''

    # Function to select a random row based on probabilities
    def select_random_row(probabilities):
        return np.random.choice(probabilities.index, p=probabilities)

    # Iterate over each row in the first DataFrame
    for index, row in df1.iterrows():
        event_name = row['name']
        probabilities = df2[event_name]
        selected_row = select_random_row(probabilities)
        
        # Fill the values in the first DataFrame
        df1.at[index, '#ref'] = df2.at[selected_row, 'window']
        df1.at[index, 'beg'] = np.random.randint(df2.at[selected_row, 'beg'], df2.at[selected_row, 'end'])
        df1.at[index, 'end'] = int(df1.at[index, 'beg']) + 1

    return(df1)


def add_elements_columns(d, df):
    '''
    Function to add the different parts of the insertions information as columns
    '''

    # create new columns
    df['Length'] = 0
    df['PolyA_Length_1'] = 0
    df['PolyA_Length_2'] = 0
    # df['Strand'] = 0
    df['RT'] = 0
    df['TSD_Length'] = 0
    df['TD_orphan_Length'] = 0
    df['VNTR_Num_Motifs'] = 0
    # df['TD_5'] = 0
    # df['TD_3'] = 0

    # iterate over rows in df
    for index, row in df.iterrows():
        name = row['name']
        for column in ['Length', 'PolyA_Length_1', 'PolyA_Length_2', 'RT', 'TSD_Length', 'TD_orphan_Length', 'VNTR_Num_Motifs']:
            key = f"{name}_{column}"
            if key in d:
                if d[key]:
                    df.loc[index, column] = d[key].pop(0)
                else:
                    df.loc[index, column] = 0
            else:
                df.loc[index, column] = 0

    # fill remaining NaN values with 0
    df.fillna(0, inplace=True)

    return(df)


def update_df(row,proportion_TD5_L1,random_numbers_dict,proportion_TD3_SVA):
    ''' 
    Function to randomly select from dictionary values and update dataframe
    '''
    if row['name'] != 'partnered_L1' and row['name'] != 'partnered_SVA':
        return pd.Series([0, 0], index=['TD_5', 'TD_3'])
    elif row['name'] == 'partnered_L1':
        random_val = random.uniform(0, 1)
        if random_val < proportion_TD5_L1:
            val = random_numbers_dict['partnered_L1_TD_5'].pop(0)
            return pd.Series([val, 0], index=['TD_5', 'TD_3'])
        else:
            val = random_numbers_dict['partnered_L1_TD_3'].pop(0)
            return pd.Series([0, val], index=['TD_5', 'TD_3'])
    elif row['name'] == 'partnered_SVA':
        random_val = random.uniform(0, 1)
        if random_val < proportion_TD3_SVA:
            val = random_numbers_dict['partnered_SVA_TD_3'].pop(0)
            return pd.Series([0, val], index=['TD_5', 'TD_3'])
        else:
            val = random_numbers_dict['partnered_SVA_TD_5'].pop(0)
            return pd.Series([val, 0], index=['TD_5', 'TD_3'])


def update_SRC_in_ref_genome(row):
    '''
    Function to update the source gene column
    '''
    if row['SRC_identifier'] == 'NA':
        row['SRC_in_ref_genome'] = 'NA'
    return row


def consensus_seqs(file_path):
    ''' 
    Function to read from a fasta file the consensus sequences

    Input: path of the fasta file
    Output: dictionary containg, for example, key: Alu_Seq - Value: the consensus sequence
    '''
    sequences = {"Alu_Seq": "", "L1_Seq": "", "SVA_Seq": ""}

    with open(file_path, "r") as file:
        lines = file.readlines()
        for i in range(len(lines)):
            if "consensus|Alu" in lines[i] and i+1 < len(lines):
                sequences["Alu_Seq"] = lines[i+1].strip()
            elif "consensus|L1" in lines[i] and i+1 < len(lines):
                sequences["L1_Seq"] = lines[i+1].strip()
            elif "consensus|SVA" in lines[i] and i+1 < len(lines):
                sequences["SVA_Seq"] = lines[i+1].strip()
            elif "NC_012920.1 Homo sapiens mitochondrion" in lines[i] and i+1 < len(lines):
                sequences["NUMT_Seq"] = lines[i+1].strip()
    return sequences



# FUNCTIONS FOR INSERTIONS

# SOLO ALU & SOLO L1
def solo_insertions(row,consensus_seqs_dict,reference_fasta):
    ''' 
    Function to generate solo_Alu and solo_L1 sequences
    '''
    RT_value = row['RT_Length']
    if RT_value == 'NA':
        return 0 
    RT_value = int(RT_value)

    if pd.isna(row['PolyA_Length_1']) or row['PolyA_Length_1'] == 'NA':
        polyA_length = 0 
    else:
        polyA_length = int(row['PolyA_Length_1'])
    
    if row['name'] == 'solo_Alu':
        seq = consensus_seqs_dict['Alu_Seq']
    elif row['name'] == 'solo_L1':
        seq = consensus_seqs_dict['L1_Seq']
    else:
        return 0
    
    # TSD
    beg = int(row['beg'])
    tsd_len = int(row['TSD_Length'])
    start_tsd = beg - tsd_len
    # Fetch a sequence using pysam
    with pysam.FastaFile(reference_fasta) as fasta_file:
        TSD = fasta_file.fetch(row['#ref'], start_tsd, beg)

    insertion = seq[-RT_value:] + 'A' * polyA_length + TSD
    return insertion


# SOLO SVA
def solo_SVA_insertions(row, dna_dict,reference_fasta):
    ''' 
    Function to generate solo_SVA sequences
    '''
    dna_seq = dna_dict['SVA_Seq']
    rt_length = int(row['RT_Length'])

    if pd.isna(row['PolyA_Length_1']) or row['PolyA_Length_1'] == 'NA':
        poly_a_length = 0 
    else:
        poly_a_length = int(row['PolyA_Length_1'])

    if pd.isna(row['SVA_Hexamer_Length']) or row['SVA_Hexamer_Length'] == 'NA':
        sva_hexamer_length = 0 
    else:
        sva_hexamer_length = int(row['SVA_Hexamer_Length'])

    total_length = int(row['Total_Length'])
    
    # Take from the DNA sequence starting from the end, as much positions as RT_Length
    truncated_seq = dna_seq[-rt_length:]
    
    # Add PolyA_Length_1 'A's to the end of the sequence
    poly_a_tail = 'A' * poly_a_length
    seq_with_poly_a = truncated_seq + poly_a_tail
    
    # Check if Total_Length is larger or equal to 1627, add SVA_Hexamer_Length CCCTCT repeats to the beginning
    if total_length >= 1565:
        sva_hexamer = 'CCCTCT' * (sva_hexamer_length // 6) + 'CCCTCT'[:sva_hexamer_length % 6]
        seq_with_sva_hexamer = sva_hexamer + seq_with_poly_a
    else:
        seq_with_sva_hexamer = seq_with_poly_a
    
    # TSD
    beg = int(row['beg'])
    tsd_len = int(row['TSD_Length'])
    start_tsd = beg - tsd_len
    # Fetch a sequence using pysam
    with pysam.FastaFile(reference_fasta) as fasta_file:
        TSD = fasta_file.fetch(row['#ref'], start_tsd, beg)

    seq_final = seq_with_sva_hexamer + TSD
    
    return seq_final



def partnered_SVA_L1_insertions(row,consensus_seqs_dict,reference_fasta):
    ''' 
    Function to generate partnered_SVA and partnered_L1 sequences
    '''
    RT_value = row['RT_Length']
    if RT_value == 'NA':
        return 0 
    RT_value = int(RT_value)

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

    # RT
    if row['name'] == 'partnered_L1':
        seq = consensus_seqs_dict['L1_Seq']
    elif row['name'] == 'partnered_SVA':
        seq = consensus_seqs_dict['SVA_Seq']
    else:
        return 0
    seq_RT = seq[-RT_value:]

    # Transduction
    # Fetch a sequence using pysam
    with pysam.FastaFile(reference_fasta) as fasta_file:
        transduction = fasta_file.fetch(row['SRC_ref'], row['TD_beg'], row['TD_end'])

    # Add SVA hexamer
    if pd.isna(row['SVA_Hexamer_Length']) or row['SVA_Hexamer_Length'] == 'NA':
        sva_hexamer_length = 0 
    else:
        sva_hexamer_length = int(row['SVA_Hexamer_Length'])

    if RT_value == 1565:
        sva_hexamer = 'CCCTCT' * (sva_hexamer_length // 6) + 'CCCTCT'[:sva_hexamer_length % 6]
    else:
        sva_hexamer = ''

    # TSD
    beg = int(row['beg'])
    tsd_len = int(row['TSD_Length'])
    start_tsd = beg - tsd_len
    # Fetch a sequence using pysam
    with pysam.FastaFile(reference_fasta) as fasta_file:
        TSD = fasta_file.fetch(row['#ref'], start_tsd, beg)
    

    return seq_RT, len(seq_RT), polyA1, len(polyA1), transduction, len(transduction), polyA2, len(polyA2), sva_hexamer, len(sva_hexamer), TSD, reference_fasta


# DUPLICATIONS
def duplication_insertions(row,reference_fasta):
    ''' 
    Function to generate the duplicated sequences
    '''
    start = int(row['beg'])
    length = int(row['Total_Length'])
    end = start + length
    with pysam.FastaFile(reference_fasta) as fasta_file:
        insertion = fasta_file.fetch(row['#ref'], start, end)
    return insertion


# INVERTED DUPLICATIONS
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



# NUMT
def numt_insertions(row,consensus_seqs_dict):
    ''' 
    Function to generate the mitochondrial insertion sequences
    '''
    # Retrieve the sequence from the dictionary
    seq = consensus_seqs_dict['NUMT_Seq']
    
    # Transform Total_Length to integer
    length = int(row['Total_Length'])
    
    # Get a random starting position from the sequence
    start_pos = random.randint(0, len(seq) - 1) #if seq else 0
    
    # Take 'length' number of positions from 'seq', wrapping around if necessary
    result = ''
    for offset in range(length):
        result += seq[(start_pos + offset) % len(seq)]
    
    return result



# ORPHAN
def orphan_insertions(row,reference_fasta):
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

    seq = transduction + polyA1 + TSD
    return seq


# VNTR
def vntr_insertions(row, separated_motifs):
    ''' 
    Function to generate VNTR sequences
    '''
    # Create a copy of the list to avoid altering the original one
    motif_list = separated_motifs.copy()
    
    # Extract the necessary values from the DataFrame row
    total_length = int(row['Total_Length'])
    number_motifs = int(row['VNTR_Num_Motifs'])
    
    selected_motifs = []
    total_selected_length = 0
    
    # Select motifs until we have the desired number of motifs and the total length is not exceeded
    while len(selected_motifs) < number_motifs:
        # Randomly select a motif
        motif = random.choice(motif_list)
        motif_length = len(motif)
        
        # Check if adding this motif would exceed the total length
        if total_selected_length + motif_length <= total_length:
            selected_motifs.append(motif)
            total_selected_length += motif_length
            motif_list.remove(motif)  # Remove the motif from the list
    
    if number_motifs == 1:
        # If only one motif, repeat until we reach total_length
        single_motif = selected_motifs[0]
        single_motif_length = len(single_motif)
        repetitions = total_length // single_motif_length
        remainder = total_length % single_motif_length
        
        sequence = single_motif * repetitions
        if remainder > 0:
            sequence += single_motif[:remainder]  # Add partial motif
            
    else:
        # If more than one motif, divide total length into non-equal parts
        lengths = random.sample(range(1, total_length // number_motifs + 1), number_motifs - 1)
        lengths.append(total_length - sum(lengths))  # Make sure total length is preserved
        
        # We ensure the lengths are varied and valid
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
    



# GLOBAL FUNCTION
def process_df(row, consensus_seqs_dict, separated_motifs,reference_fasta):
    '''
    Global function that for each row of the df generates the determied sequence
    '''
    # SOLO insertions
    if row['name'] in ['solo_Alu', 'solo_L1']:
        return solo_insertions(row,consensus_seqs_dict,reference_fasta)
    elif row['name'] in ['solo_SVA']:
        return solo_SVA_insertions(row, consensus_seqs_dict,reference_fasta)

    # Partnered SVA
    elif row['name'] in ['partnered_SVA']:
        results = partnered_SVA_L1_insertions(row,consensus_seqs_dict,reference_fasta)
        # TD_5_Length != 0 means it has 5' transduction
        if row['TD_5_Length'] != 'NA':
            # Transduction + Hexamer + RT + polyA1
            return results[4] + results[8] + results[0] + results[2] + results[10]
        # If not it has 3' transduction
        else:
            # HEXAMER + RT + PolyA1 + Transcution + PolyA2
            return results[8] + results[0] + results[2] + results[4] + results[6] + results[10]
    
    # Partnered L1
    elif row['name'] in ['partnered_L1']:
        results = partnered_SVA_L1_insertions(row,consensus_seqs_dict,reference_fasta)
        # TD_5_Length != 0 means it has 5' transduction
        if row['TD_5_Length'] != 'NA':
            # Tranduction + RT + PolyA1
            return results[4] + results[0] + results[2] + results[10]
        # If not it has 3' transduction
        else:
            # RT + PolyA1 + Transduction + PolyA2
            return results[0] + results[2] + results[4] + results[6] + results[10]
    
    # DUPLICATIONS
    elif row['name'] == 'DUP':
        return duplication_insertions(row,reference_fasta)     
    
    # INVERSE DUPLICATIONS  
    elif row['name'] == 'INV_DUP':
        seq = duplication_insertions(row) 
        return reverse_complementary(seq,reference_fasta)
    
    # NUMT
    elif row['name'] == 'NUMT':
        return numt_insertions(row,consensus_seqs_dict)
    
    # ORPHAN
    elif row['name'] in ['orphan']:
        return orphan_insertions(row,reference_fasta)   

    # VNTR
    elif row['name'] in ['VNTR']:
        return vntr_insertions(row, separated_motifs) 
    else:
        return 0


def RC_insertion(df):
    '''
    Function to create the reverse complementary of those new generated sequences in the negative strand
    '''
    df['sequence_ins'] = df.apply(
        lambda row: reverse_complementary(row['sequence_ins']) if row['Strand'] == '-' else row['sequence_ins'],
        axis=1
    )
    return df


def main(vcf_path_insertions, chromosome_length, bin_size):
    print(f'VCF file with insertions: {vcf_path_insertions}')
    print(f'File with chromosomes length: {chromosome_length}')

    # From the VCF create directly a BED object
    result_df = read_vcf_file_BED(vcf_path_insertions)

    # Change the names to the correct ones:
    result_BED_table = result_df.rename(columns={'Type_SV': 'name', 'Chromosome':'#ref', 'Start_position': 'beg', 'End_position':'end', 'Family':'SubType'})

    # Replace all NaN values with 'NA'
    result_BED_table.fillna('NA', inplace=True)

    # Join name and SubType by '_'
    result_BED_table['name'] = result_BED_table['name'] + '_' + result_BED_table['SubType'].apply(lambda x: x if x in ['SVA', 'Alu', 'L1'] else '')
    # Delete the SubType column
    result_BED_table = result_BED_table.drop('SubType', axis=1)
    # Remove '_' from elements that end with '_' 
    result_BED_table['name'] = result_BED_table['name'].apply(lambda x: x[:-1] if x.endswith('_') else x)


    # Applying the filter to get only the ones that have 1 transduction
    TD_filter(result_BED_table)

    # From the PolyA_Length lets get the 1st value in one column and the second value in another column
    result_BED_table['PolyA_Length_1'] = result_BED_table['PolyA_Length'].apply(lambda x: x.split(',')[0] if x != 'NA' else x)
    result_BED_table['PolyA_Length_2'] = result_BED_table['PolyA_Length'].apply(lambda x: x.split(',')[1] if x != 'NA' and len(x.split(',')) > 1 else 'NA')

    # Delete the original PolyA_Length column
    result_BED_table.drop('PolyA_Length', axis=1, inplace=True)

    # convert 'PolyA_Length' columns to numeric, converting 'NA' to NaN
    result_BED_table['PolyA_Length_1'] = pd.to_numeric(result_BED_table['PolyA_Length_1'], errors='coerce')
    result_BED_table['PolyA_Length_2'] = pd.to_numeric(result_BED_table['PolyA_Length_2'], errors='coerce')

    # Manually calculate RT lengh
    result_BED_table.drop('RT_Length', axis=1, inplace=True)

    # Transform TSD_Length & SVA_Hexamer column to numerical
    result_BED_table['TSD_Length'] = pd.to_numeric(result_BED_table['TSD_Length'], errors='coerce')
    result_BED_table['SVA_Hexamer'] = pd.to_numeric(result_BED_table['SVA_Hexamer'], errors='coerce')

    # create the new 'RT' column and add the values
    result_BED_table['RT'] = result_BED_table.apply(lambda row: row['Length'] - row['TSD_Length'] - row['PolyA_Length_1'] - row['PolyA_Length_2'] - (row['SVA_Hexamer'] if pd.notnull(row['SVA_Hexamer']) else 0) if pd.notnull(row['PolyA_Length_2']) else row['Length'] - row['TSD_Length'] - row['PolyA_Length_1'] - (row['SVA_Hexamer'] if pd.notnull(row['SVA_Hexamer']) else 0) if pd.notnull(row['PolyA_Length_1']) else None, axis=1)

    # Create dictionary
    dict_mutations = dict_mutations_list(result_BED_table)

    # Now let's remove the non desired keys:
    keys_remove = ['chimera_PolyA_Length', 'chimera_TSD_Length', 'chimera_Length',
                    'partnered_PolyA_Length', 'partnered_RT', 'partnered_TSD_Length',
                    'partnered_Length', 'partnered_Alu_PolyA_Length', 'partnered_Alu_RT',
                    'partnered_Alu_TSD_Length', 'partnered_Alu_Length','chimera_RT','partnered_RT','partnered_PolyA_Length','partnered_RT','partnered_TSD_Length',
                    'partnered_Length', 'VNTR_PolyA_Length', 'VNTR_RT', 'DUP_PolyA_Length','DUP_RT', 'INV_DUP_PolyA_Length','INV_DUP_RT','COMPLEX_DUP_PolyA_Length',
                    'COMPLEX_DUP_RT','NUMT_PolyA_Length','NUMT_RT','DUP_INTERSPERSED_Alu_PolyA_Length','DUP_INTERSPERSED_Alu_RT','DUP_INTERSPERSED_L1_PolyA_Length',
                    'DUP_INTERSPERSED_L1_RT','DUP_INTERSPERSED_SVA_PolyA_Length', 'DUP_INTERSPERSED_SVA_RT','DUP_INTERSPERSED_PolyA_Length','DUP_INTERSPERSED_RT',
                    'COMPLEX_DUP_Length','chimera_TD_5','chimera_TD_3','partnered_TD_5','partnered_TD_3','partnered_Alu_TD_5','partnered_Alu_TD_3',
                    'solo_SVA_PolyA_Length_2','solo_Alu_PolyA_Length_2','partnered_Alu_PolyA_Length_1','partnered_Alu_PolyA_Length_2','orphan_PolyA_Length_2',
                    'solo_L1_PolyA_Length_2','PSD_PolyA_Length_2','chimera_PolyA_Length_2','chimera_PolyA_Length_1','partnered_PolyA_Length_1','partnered_PolyA_Length_2',
                    'VNTR_PolyA_Length_1','VNTR_PolyA_Length_2','DUP_PolyA_Length_1','DUP_PolyA_Length_2','INV_DUP_PolyA_Length_1','INV_DUP_PolyA_Length_2',
                    'COMPLEX_DUP_PolyA_Length_1','COMPLEX_DUP_PolyA_Length_2','NUMT_PolyA_Length_1','NUMT_PolyA_Length_2','DUP_INTERSPERSED_Alu_PolyA_Length_1',
                    'DUP_INTERSPERSED_Alu_PolyA_Length_2','DUP_INTERSPERSED_L1_PolyA_Length_1','DUP_INTERSPERSED_L1_PolyA_Length_2','DUP_INTERSPERSED_SVA_PolyA_Length_1',
                    'DUP_INTERSPERSED_SVA_PolyA_Length_2','DUP_INTERSPERSED_PolyA_Length_1','DUP_INTERSPERSED_PolyA_Length_2','VNTR_TSD_Length','DUP_TSD_Length',
                    'INV_DUP_TSD_Length','COMPLEX_DUP_TSD_Length','NUMT_TSD_Length','DUP_INTERSPERSED_Alu_TSD_Length','DUP_INTERSPERSED_L1_TSD_Length',
                    'DUP_INTERSPERSED_SVA_TSD_Length','DUP_INTERSPERSED_TSD_Length','solo_Alu_SVA_Hexamer','partnered_Alu_SVA_Hexamer','partnered_L1_SVA_Hexamer',
                    'orphan_SVA_Hexamer','solo_L1_SVA_Hexamer','PSD_SVA_Hexamer','chimera_SVA_Hexamer','partnered_SVA_Hexamer','VNTR_SVA_Hexamer','DUP_SVA_Hexamer',
                    'INV_DUP_SVA_Hexamer','COMPLEX_DUP_SVA_Hexamer','NUMT_SVA_Hexamer','DUP_INTERSPERSED_Alu_SVA_Hexamer','DUP_INTERSPERSED_L1_SVA_Hexamer',
                    'DUP_INTERSPERSED_SVA_SVA_Hexamer','DUP_INTERSPERSED_SVA_Hexamer','DUP_INTERSPERSED_Alu_Length','DUP_INTERSPERSED_L1_Length','DUP_INTERSPERSED_SVA_Length','DUP_INTERSPERSED_Length',
                    'PSD_PolyA_Length_1' , 'PSD_RT' , 'PSD_TSD_Length' , 'PSD_Length']

    for key in keys_remove:
        dict_mutations.pop(key, None)

    # Remove outliers
    keys_filter = ['INV_DUP_Length','DUP_Length','VNTR_Length']
    filter_sd(dict_mutations, keys_filter)

    # Specify the keys for which you want to delete nan values
    keys_remove_nan = ['partnered_SVA_PolyA_Length_2', 'partnered_L1_PolyA_Length_2','solo_SVA_SVA_Hexamer','partnered_SVA_SVA_Hexamer', 'TD_orphan_Length', 'VNTR_Num_Motifs', 'VNTR_Motifs']

    # Loop through the specified keys and remove nan values
    for key in keys_remove_nan:
        if key in dict_mutations:
            dict_mutations[key] = remove_nan(dict_mutations[key])

    # Important for VNTR
    # Here we extract from the dictionary the motifs of VNTR (the letters) and we keep them in a list, and remove them from the dictionary
    # Extract the values from the key 'VNTR_VNTR_Motifs'
    if 'VNTR_VNTR_Motifs' in dict_mutations:
        # Get the list of values associated with the key
        motifs = dict_mutations['VNTR_VNTR_Motifs']
        
        # Create a new list to hold the processed values
        separated_motifs = []

        # Iterate through each item in the list
        for motif in motifs:
            # Split the string by ',' and extend the separated_motifs list
            separated_motifs.extend(motif.split(','))

        # Remove the key from the dictionary
        del dict_mutations['VNTR_VNTR_Motifs']

        # Now `separated_motifs` contains the separated values



    # Remove NA's and transform to integrers all values
    for key, value_list in dict_mutations.items():
        new_value_list = []
        for val in value_list:
            if val == 'NA' or isinstance(val, float) and math.isnan(val):
                new_value_list.append(val)
            else:
                new_value_list.append(int(val))
        dict_mutations[key] = new_value_list


    dict_mutations = {key: value for key, value in dict_mutations.items() if not isinstance(value, float) or not math.isnan(value)}



    ##########################################
    # TABLE INSERTIONS' FEATURES
    #########################################
    # Define rows and columns
    rows = ['VNTR', 'solo_Alu', 'DUP', 'solo_L1', 'partnered_L1', 'solo_SVA', 'partnered_SVA', 'orphan', 'NUMT', 'INV_DUP']
    columns = ['PolyA_Length_1', 'RT', 'TSD_Length', 'Length', 'SVA_Hexamer', 'PolyA_Length_2', 'TD_3', 'TD_5', 'VNTR_Num_Motifs', 'TD_orphan_Length']

    # Create a DataFrame initialized with None
    insertion_features_df = pd.DataFrame(index=rows, columns=columns)

    # Populate the DataFrame with the values from the dictionary
    for key, values in dict_mutations.items():
        # Split the key to get the type and the measurement
        key_parts = key.split('_')

        # Determine the row name based on the key
        if key.startswith('solo'):
            row_name = 'solo_' + key_parts[1]
        elif key.startswith('partnered'):
            row_name = 'partnered_' + key_parts[1]
        elif key.startswith('orphan_'):
            row_name = 'orphan'
        elif key.startswith('INV_'):
            row_name = 'INV_DUP'
        else:
            row_name = key_parts[0]  # for 'VNTR' or 'DUP'
        
        # Determine the column name based on the key parts
        if len(key_parts) > 2:
            column_name = "_".join(key_parts[2:])  # Get everything after the second part as column name
        else:
            column_name = key_parts[1]  # In case the key only has two parts, use the second part for column name

        # Fill the DataFrame if row and column names match
        if row_name in insertion_features_df.index and column_name in insertion_features_df.columns:
            insertion_features_df.at[row_name, column_name] = values

    # Special case for VNTR_Num_Motifs 
    if 'VNTR_VNTR_Num_Motifs' in dict_mutations:
        insertion_features_df.at['VNTR', 'VNTR_Num_Motifs'] = dict_mutations['VNTR_VNTR_Num_Motifs']

    # Now handle orphan data separately
    for key, values in dict_mutations.items():
        if key.startswith('orphan_'):
            # Determine the column name based on the key parts
            column_name = "_".join(key.split('_')[1:])  # Skip the 'orphan_' prefix
            if column_name in insertion_features_df.columns:
                insertion_features_df.at['orphan', column_name] = values

    # Function to convert list cells to strings without brackets
    def convert_list_to_string(value):
        if isinstance(value, list):
            return ','.join(map(str, value))  # Join numbers by commas
        return value  # Leave other types unchanged

    insertion_features_df = insertion_features_df.applymap(convert_list_to_string)

    insertion_features_df.to_csv('Insertion_Features.tsv', sep='\t', index=True)

    ###################################################
    # GENOME WIDE DISTRIBUTION
    ###################################################

    # Distribution of events genome-wide
    # Dictionary containing references as keys and their lengths as values:
    chr_length = formats.chrom_lengths_index(chromosome_length)

    # Bin size 
    binSize = bin_size # this is 1MB
    bins = gRanges.makeGenomicBins(chr_length, binSize, ['chr1','chr2','chr3','chr4','chr5','chr6','chr7','chr8','chr9','chr10','chr11','chr12','chr13','chr14','chr15','chr16','chr17','chr18','chr19','chr20','chr21','chr22'])[::-1] # defino tamaño ventana y con las longitudes y tamanño ventana crea las ventanas 

    # Create the table of insertions classified in windows
    res_table = mut_bins(bins, result_BED_table)

    prob_INS = calculate_normalized_values(res_table)

    # Replace all NaN values with 0.00
    prob_INS.fillna(0.00, inplace=True)

    # Now let's delete the columns of the variants we are not interested:
    columns_to_delete = ['partnered_Alu', 'chimera', 'partnered', 'COMPLEX_DUP', 'PSD']
    # Delete the columns
    prob_INS = prob_INS.drop(columns=columns_to_delete)

    # Table of probabilities normalized by column:
    prob_INS_col_normalized = normalize_columns(prob_INS)

    prob_INS_col_normalized.to_csv('Genome_Wide_Distribution.tsv', sep='\t', index=False)

    ################################################
    # PROBABILITIES
    ################################################

    # Get how many times each event is in the VCF table
    name_distribution = result_BED_table['name'].value_counts()

    # Drop the names we are not interested in
    name_distribution = name_distribution.drop(['partnered_Alu', 'chimera', 'partnered', 'COMPLEX_DUP','DUP_INTERSPERSED','DUP_INTERSPERSED_Alu','DUP_INTERSPERSED_L1','DUP_INTERSPERSED_SVA', 'PSD'])

    # Transform it to a df
    name_distribution_df = pd.DataFrame({'Event': name_distribution.index, 'number': name_distribution.values})

    # Calculate probabilities
    total = name_distribution_df['number'].sum()
    name_distribution_df['Probability'] = name_distribution_df['number'] / total
    name_distribution_df = name_distribution_df.drop(columns=['number'])

    name_distribution_df.to_csv('Probabilities.tsv', sep='\t', index=False)



if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Process insertions from VCF to BED format.')
    parser.add_argument('vcf_path_insertions', type=str, help='Path to the VCF file containing insertions.')
    parser.add_argument('chromosome_length', type=str, help='Path to the chromosome length file.')
    parser.add_argument('--bin_size', type=int, default=1000000, help='Size of genomic bins (default: 1000000).')
    args = parser.parse_args()
    main(args.vcf_path_insertions, args.chromosome_length, args.bin_size)