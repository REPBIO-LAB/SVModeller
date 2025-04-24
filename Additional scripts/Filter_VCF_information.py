# SVModeller 

# Additional script to remove information from VCF file keeping just the event length

# Input:
# - VCF input file
# - VCF output file

# Output:
# - Modified VCF file without information field, just keeping the event length

import argparse

def filter_vcf_info(input_vcf, output_vcf):
    with open(input_vcf, 'r') as infile, open(output_vcf, 'w') as outfile:
        for line in infile:
            if line.startswith('#'):
                outfile.write(line)
                continue

            columns = line.strip().split('\t')
            info_field = columns[7]

            # Filter INFO fields to keep only INS_LEN or DEL_LEN
            info_parts = info_field.split(';')
            filtered_info = [item for item in info_parts if item.startswith('INS_LEN=') or item.startswith('DEL_LEN=')]

            columns[7] = ';'.join(filtered_info) if filtered_info else '.'
            outfile.write('\t'.join(columns) + '\n')
        
        print('VCF generated successfully')

def main(input_vcf, output_vcf):
    # Print the paths of the input files
    print(f'VCF file path: {input_vcf}')
    print(f'VCF output path: {output_vcf}')
    
    # Ensure .vcf extension
    if not output_vcf.endswith('.vcf'):
        output_vcf += '.vcf'
    
    filter_vcf_info(input_vcf, output_vcf)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Process VCF file to remove all information except from the length of the events')
    parser.add_argument('input_vcf', type=str, help='Path to the VCF.')
    parser.add_argument('output_vcf', type=str, help='Resulting VCF without information name.')

    args = parser.parse_args()
    main(args.input_vcf, args.output_vcf)
