# SVModeler

# Module 1.2


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


def generate_dict_from_table(data):
    result_dict = {}
    for index, row in data.iterrows():
        event = row['Event']
        for column in data.columns[1:]:  # Skip the 'Event' column
            value = row[column]
            if pd.notna(value):  # Check if value is not NaN
                key = f"{event}_{column}"
                result_dict[key] = list(map(int, value.split(',')))
    return result_dict





def main(genome_wide_df, insertion_features_df, probabilities_numbers_df, table_source_L1, table_source_SVA, chromosome_length, consensus_path, reference_fasta, num_events):
    print(f'File with genome-wide distribution: {genome_wide_df}')
    print(f'File with insertion features: {insertion_features_df}')
    print(f'File with probabilities or defined number of events: {probabilities_numbers_df}')
    print(f'File with source loci for SVA transductions: {table_source_SVA}')
    print(f'File with source loci for LINE-1 transductions: {table_source_L1}')
    print(f'File with chromosomes length: {chromosome_length}')
    print(f'File with consensus sequences: {consensus_path}')
    print(f'File with reference genome: {reference_fasta}')
    print(f'Defined number of events: {num_events}')


    table = pd.read_csv(probabilities_numbers_df, sep='\t')

    if 'Probability' in table.columns:
        sampled_names = np.random.choice(table['Event'], size=num_events, p=table['Probability'])
        df_insertions = pd.DataFrame({'name': sampled_names})
    elif 'Number' in table.columns:
        rows = []
        for i, row in table.iterrows():
            event = row['Event']
            count = row['Number']
            rows.extend([event] * count)
        df_insertions = pd.DataFrame({'name': rows})
    else:
        raise ValueError("The second column must be either 'Probability' or 'Number'.")


    table_insertion_features = pd.read_csv(insertion_features_df, sep='\t')
    table_insertion_features.columns.values[0] = 'Event'

    dict_insertion_features = generate_dict_from_table(table_insertion_features)



    # As to run the package for the distributions I need np arrays, let's create a dictionary as the previous one but with np arrays.
    # empty dict
    array_dict = {}

    for key in dict_insertion_features:
        # create a numpy array for each key
        array = np.array(dict_insertion_features[key])
        array_dict[key] = array


    # From the dictionary, of the values of the VCF, generate the distributions, and from these distributions generate a dictionary with determined random numbers based on the distributions.
    number_random_events=num_events*3
    random_numbers_dict = generate_random_numbers(array_dict, number_random_events)

    # As there are some negative values from the distribution, let's remove them
    # Iterate through the dictionary and remove negative values
    for key, values in random_numbers_dict.items():
        random_numbers_dict[key] = remove_negative_values(values)

    # Now let's remove those elements from the keys solo_Alu_RT, solo_L1_RT the values that are higher than the length of the consensus sequence of the element.
    # Because if not, random numbers at the tail of the graph are also generated.
    # dictionary with length of each consensus sequence
    length_consensus_seq = {
        'solo_Alu_RT': 281,
        'partnered_Alu_RT': 281,
        'solo_L1_RT': 6023,
        'partnered_L1_RT': 6023,
        'solo_SVA_RT': 1565,
        'partnered_SVA_RT': 1565,
    }

    for key, value in length_consensus_seq.items():
        if key in random_numbers_dict:
            random_numbers_dict[key] = [x for x in random_numbers_dict[key] if x <= value]



    # Now I take the df of the names in each row AND I add a column beg end and ref based on the probabilities table
    df_insertions2 = add_beg_end_columns(df_insertions,genome_wide_df)

    # Reorder the columns so I create a BED file order of the columns
    df_insertions2 = df_insertions2[['#ref', 'beg', 'end', 'name']]

    df_insertions3 = add_elements_columns(random_numbers_dict, df_insertions2)

    # ### Proportions of Transductions
    # Get how many partnered SVA and L1 have transductions


    check_TD = result_BED_table
    check_TD['TD_5'] = pd.to_numeric(check_TD['TD_5'], errors='coerce')
    check_TD['TD_3'] = pd.to_numeric(check_TD['TD_3'], errors='coerce')
    # count number of 'partnered_SVA' in the 'column' column
    partnered_SVA_count = len(check_TD[check_TD['name'] == 'partnered_SVA'])

    # count number of 'partnered_SVA' with non-zero values in 'TD_5' column
    partnered_SVA_TD5_nonzero_count = len(check_TD[(check_TD['name'] == 'partnered_SVA') & (check_TD['TD_5'] > 0)])

    # count number of 'partnered_SVA' with non-zero values in 'TD_3' column
    partnered_SVA_TD3_nonzero_count = len(check_TD[(check_TD['name'] == 'partnered_SVA') & (check_TD['TD_3'] > 0)])

    # count number of 'partnered_L1' in the 'column' column
    partnered_L1_count = len(check_TD[check_TD['name'] == 'partnered_L1'])

    # count number of 'partnered_L1' with non-zero values in 'TD_5' column
    partnered_L1_TD5_nonzero_count = len(check_TD[(check_TD['name'] == 'partnered_L1') & (check_TD['TD_5'] > 0)])

    # count number of 'partnered_L1' with non-zero values in 'TD_3' column
    partnered_L1_TD3_nonzero_count = len(check_TD[(check_TD['name'] == 'partnered_L1') & (check_TD['TD_3'] > 0)])

    # Total values TD SVA
    total_TD_SVA = partnered_SVA_TD5_nonzero_count + partnered_SVA_TD3_nonzero_count
    # Proportions
    proportion_TD5_SVA = partnered_SVA_TD5_nonzero_count / total_TD_SVA
    proportion_TD3_SVA = partnered_SVA_TD3_nonzero_count / total_TD_SVA

    # Total values TD SVA
    total_TD_L1 = partnered_L1_TD5_nonzero_count + partnered_L1_TD3_nonzero_count
    # Proportions
    proportion_TD5_L1 = partnered_L1_TD5_nonzero_count / total_TD_L1
    proportion_TD3_L1 = partnered_L1_TD3_nonzero_count / total_TD_L1


    # Call the update_df function with extra arguments using lambda
    df_insertions3[['TD_5', 'TD_3']] = df_insertions3.apply(lambda row: update_df(row, proportion_TD5_L1, random_numbers_dict, proportion_TD3_SVA), axis=1)

    # Now let's remove the second polyA from the Transduction 5' ones
    df_insertions3.loc[df_insertions3['TD_5'] != 0, 'PolyA_Length_2'] = 0



    # Add Hexamer Length to the SVA's
    # Create a new column 'SVA_Hexamer_Length' in the dataframe
    df_insertions3['SVA_Hexamer_Length'] = 0

    for index, row in df_insertions3.iterrows():
        if row['name'] == 'solo_SVA':
            if row['RT'] >= 1565:
                df_insertions3.at[index, 'RT'] = 1565
                df_insertions3.at[index, 'SVA_Hexamer_Length'] = random_numbers_dict['solo_SVA_SVA_Hexamer'][0]
                random_numbers_dict['solo_SVA_SVA_Hexamer'] = random_numbers_dict['solo_SVA_SVA_Hexamer'][1:]
            elif row['RT'] < 1565:
                df_insertions3.at[index, 'SVA_Hexamer_Length'] = 0
        elif row['name'] == 'partnered_SVA':
            if row['RT'] >= 1565:
                df_insertions3.at[index, 'RT'] = 1565
                df_insertions3.at[index, 'SVA_Hexamer_Length'] = random_numbers_dict['partnered_SVA_SVA_Hexamer'][0]
                random_numbers_dict['partnered_SVA_SVA_Hexamer'] = random_numbers_dict['partnered_SVA_SVA_Hexamer'][1:]
            elif row['RT'] <= 1565:
                df_insertions3.at[index, 'SVA_Hexamer_Length'] = 0


    # #### Total Length column update
    # First let's change the name of the column Length to Total_Length
    df_insertions3 = df_insertions3.rename(columns={'Length': 'Total_Length'})

    # Update column
    conditions = ['solo_Alu', 'partnered_L1', 'solo_SVA', 'partnered_SVA', 'solo_L1', 'orphan', 'PSD']
    columns_to_sum = ['PolyA_Length_1', 'PolyA_Length_2', 'RT', 'TSD_Length', 'TD_5', 'TD_3','SVA_Hexamer_Length', 'TD_orphan_Length']

    for index, row in df_insertions3.iterrows():
        if row['name'] in conditions:
            df_insertions3.at[index, 'Total_Length'] = row[columns_to_sum].sum()

    # #### Add column strand
    # Add column strand of the insertion
    df_insertions3['Strand'] = [random.choice(['+', '-']) for x in range(len(df_insertions3))]

    # #### Update name of columns
    df_insertions3.rename(columns={'RT': 'RT_Length', 'TD_5': 'TD_5_Length', 'TD_3': 'TD_3_Length'}, inplace=True)


    # ## Add the source element information & calculations
    # - 1st add a row to the insertions table based on the probability of the source element to the partnered SVA and L1 elements
    # - Do calculations of the transduction length based on the string of the source element

    # Add columns from the first dataframe to the second dataframe
    df_insertions3[['SRC_identifier', 'SRC_ref', 'SRC_beg', 'SRC_end', 'SRC_cont_PCAWG', 'SRC_strand', 'SRC_in_ref_genome']] = 0
    table_source_L1 = pd.read_csv(args.table_source_L1, sep='\t')
    table_source_SVA = pd.read_csv(args.table_source_SVA, sep='\t')
    # Iterate over rows in the second dataframe
    for index, row in df_insertions3.iterrows():
        if row['name'] in ['partnered_SVA', 'partnered_L1', 'orphan']:
            # Select a row from the first dataframe based on probabilities in SRC_cont_PCAWG
            probabilities = table_source_L1['SRC_cont_PCAWG'].values
            selected_row = table_source_L1.sample(weights=probabilities).iloc[0]
            # Add the selected row to the second dataframe
            df_insertions3.loc[index, ['SRC_identifier', 'SRC_ref', 'SRC_beg', 'SRC_end', 'SRC_cont_PCAWG', 'SRC_strand', 'SRC_in_ref_genome']] = selected_row.values

    df_insertions3.fillna(0, inplace=True)


    # Add columns from the first dataframe to the second dataframe
    df_insertions3[['SRC_identifier', 'SRC_ref', 'SRC_beg', 'SRC_end', 'SRC_cont_PCAWG', 'SRC_contribution', 'SRC_strand', 'SRC_in_ref_genome']] = 0

    # Iterate over rows in the second dataframe
    for index, row in df_insertions3.iterrows():
        if row['name'] in ['partnered_L1', 'orphan']:
            # Select a row from the first dataframe based on probabilities in SRC_cont_PCAWG
            probabilities = table_source_L1['SRC_cont_PCAWG'].values
            selected_row = table_source_L1.sample(weights=probabilities).iloc[0]
            # Add the selected row to the second dataframe
            df_insertions3.loc[index, ['SRC_identifier', 'SRC_ref', 'SRC_beg', 'SRC_end', 'SRC_cont_PCAWG',  'SRC_strand', 'SRC_in_ref_genome']] = selected_row.values
        elif row['name'] in ['partnered_SVA']:
            # Select a row from the first dataframe based on probabilities in SRC_contribution
            probabilities = table_source_SVA['SRC_contribution'].values
            selected_row = table_source_SVA.sample(weights=probabilities).iloc[0]
            # Add the selected row to the second dataframe
            df_insertions3.loc[index, ['SRC_identifier', 'SRC_ref', 'SRC_beg', 'SRC_end','SRC_contribution', 'SRC_strand', 'SRC_in_ref_genome']] = selected_row.values
    df_insertions3.fillna(0, inplace=True)

    # Add new columns 'TD_beg' and 'TD_end' TD_orphan_Length
    df_insertions3['TD_beg'] = 0
    df_insertions3['TD_end'] = 0

    # Modify the values in 'TD_beg' and 'TD_end' based on conditions
    for index, row in df_insertions3.iterrows():
        if row['name'] == 'partnered_SVA' or row['name'] == 'partnered_L1':
            if row['SRC_strand'] == 'plus':
                df_insertions3.at[index, 'TD_beg'] = row['SRC_end']
                df_insertions3.at[index, 'TD_end'] = row['SRC_end'] + row['TD_5_Length'] if row['TD_5_Length'] != 0 else row['SRC_end'] + row['TD_3_Length']
            elif row['SRC_strand'] == 'minus':
                df_insertions3.at[index, 'TD_end'] = row['SRC_beg']
                df_insertions3.at[index, 'TD_beg'] = row['SRC_beg'] - row['TD_5_Length'] if row['TD_5_Length'] != 0 else row['SRC_beg'] - row['TD_3_Length']
        if row['name'] == 'orphan':
            if row['SRC_strand'] == 'plus':
                df_insertions3.at[index, 'TD_beg'] = row['SRC_end']
                df_insertions3.at[index, 'TD_end'] = row['SRC_end'] + row['TD_orphan_Length']
            elif row['SRC_strand'] == 'minus':
                df_insertions3.at[index, 'TD_end'] = row['SRC_beg']
                df_insertions3.at[index, 'TD_beg'] = row['SRC_beg'] - row['TD_orphan_Length']


    # Modify SRC_strand
    df_insertions3['SRC_strand'] = df_insertions3['SRC_strand'].replace({'minus': '-', 'plus': '+'})


    # Change all the 0s for 'NA'
    # Specify the columns where 0s should not be transformed to NA
    keep_columns = ['SRC_in_ref_genome']
    # Transform all the 0s to NA except for the specified columns
    for column in df_insertions3.columns:
        if column not in keep_columns:
            df_insertions3[column] = df_insertions3[column].apply(lambda x: 'NA' if x == 0 else x)


    # In the column SRC_in_ref_genome now let's change those 0s that are of the elements that do not have a source elemnt insertion
    # Apply the function
    df_insertions3 = df_insertions3.apply(update_SRC_in_ref_genome, axis=1)

    # Get the consensus sequences
    consensus_seqs_dict = consensus_seqs(consensus_path)
        
    # Generate column with final insertion sequence
    df_insertions3['sequence_ins'] = df_insertions3.apply(lambda row: process_df(row, consensus_seqs_dict,separated_motifs,reference_fasta), axis=1)

    # Add a new column 'event_type' with the value 'insertion'
    df_insertions3['event_type'] = 'insertion'

    # Reorder columns to make 'event_type' the 4th column
    # Create a new column order
    new_order = df_insertions3.columns.tolist()
    new_order.remove('event_type')  # Remove 'event_type' from the current position
    new_order.insert(3, 'event_type')  # Insert 'event_type' as the 4th column

    # Reassign the DataFrame to reorder the columns
    df_insertions3 = df_insertions3[new_order]

    # Write all the insertion sequence in Upper case
    df_insertions3['sequence_ins'] = df_insertions3['sequence_ins'].str.upper()

    # Apply reverse complementary
    df_insertions3 = RC_insertion(df_insertions3)

    # Save the output
    df_insertions3.to_csv('Insertions_table.tsv', sep='\t', index=False)




if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Process insertions from VCF to BED format.')
    parser.add_argument('genome_wide_df', type=str, help='Path to the TSV file containing genome-wide distribution of events.')
    parser.add_argument('insertion_features_df', type=str, help='Path to the TSV file containing insertion features of events.')
    parser.add_argument('probabilities_numbers_df', type=str, help='Path to the TSV file probabilities or defined number of events.')
    parser.add_argument('table_source_L1', type=str, help='Path to the TSV file containing loci for LINE-1 transductions.')
    parser.add_argument('table_source_SVA', type=str, help='Path to the TSV file containing loci for SVA transductions.')
    parser.add_argument('chromosome_length', type=str, help='Path to the chromosome length file.')
    parser.add_argument('consensus_path', type=str, help='Path to file with consensus sequences.')
    parser.add_argument('reference_fasta', type=str, help='Path to file with reference genome.')
    parser.add_argument('--num_events', type=int, default=100, help='Number of events to sample (optional, just in case of providing probabilities).')
    args = parser.parse_args()
    main(args.genome_wide_df, args.insertion_features_df, args.probabilities_numbers_df, args.table_source_L1, args.table_source_SVA, args.chromosome_length, args.consensus_path, args.reference_fasta, args.num_events)
