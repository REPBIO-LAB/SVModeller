# SVModeller - Module 5

# Generate BAM file with reads at different coverage and allele frequency levels

# Input:
# - Reference genome (.fasta)
# - Modified reference genome (.fasta)
# - Method to generate reads (choices='quality_score', 'error_model', 'training')
# - Method file (string, '.model' or '.fastq' i.e. ERRHMM-ONT-HQ.model)
# - Allele frequency (decimal number between 0-1. i.e. 0.25)
# - Coverage (integer number: 30 for 30x)
# - Output directory (file name)
# - Technology of generated reads (ONT, PB, HiFi)

# Output:
# - Alignment (combined_final_alignment.bam)
# - Modified genome reads (.fastq)
# - Reference genome reads (.fastq)

# Developers
# SVModeller has been developed by Ismael Vera-Munoz (orcid.org/0009-0009-2860-378X) at the Repetitive DNA Biology (REPBIO) Lab at the Centre for Genomic Regulation (CRG) (Barcelona 2024-2026)

# License
# SVModeller is distributed under the AGPL-3.0.

import subprocess
import os
import argparse
import glob
import warnings

# Function to run PBSIM to generate synthetic reads for reference and modified genomes
def run_pbsim(genome, method_file, method, depth, output_dir, output_reference):
    if depth == 0:
        print("Depth is 0. Skipping PBSIM execution.")
        return
    
    method_file = os.path.abspath(method_file)
    output_prefix = os.path.join(output_dir, output_reference)

    if method == 'quality_score':
        command = f"pbsim --strategy wgs --method qshmm --qshmm {method_file} --depth {depth} --genome {genome} --prefix {output_prefix}"
    elif method == 'error_model':
        command = f"pbsim --strategy wgs --method errhmm --errhmm {method_file} --depth {depth} --genome {genome} --prefix {output_prefix}"
    elif method == 'training':
        command = f"pbsim --strategy wgs --method sample --sample {method_file} --depth {depth} --genome {genome} --prefix {output_prefix}"
    else:
        raise ValueError(f"Unknown method: {method}")

    print(f"Running PBSIM with command: {command}")
    subprocess.run(command, shell=True, check=True)

    return output_prefix

# Function to align reads using Minimap2
def run_minimap2(reference_file, fastq_file_1, fastq_file_2, output_bam, technology, threads):
    fastqs = fastq_file_1
    if fastq_file_2:
        fastqs += f" {fastq_file_2}"

    if technology == 'ONT':
        command = f"minimap2 -ax map-ont {reference_file} {fastqs} -t {threads}"
    elif technology == 'PB':
        command = f"minimap2 -ax map-pb {reference_file} {fastqs} -t {threads}"
    elif technology == 'HiFi':
        command = f"minimap2 -ax map-hifi {reference_file} {fastqs} -t {threads}"
    else:
        raise ValueError(f"Unknown technology: {technology}")

    command += f" | samtools view -bS -o {output_bam} -@ {threads}"
    subprocess.run(command, shell=True, check=True)

# Function to sort BAM file
def sort_bam(bam_file, threads):
    sorted_bam_file = bam_file.replace('.bam', '.sorted.bam')
    command = f"samtools sort {bam_file} -o {sorted_bam_file} -@ {threads}"
    print(f"Sorting BAM file with command: {command}")
    subprocess.run(command, shell=True, check=True)
    return sorted_bam_file

# Function to index BAM file
def index_bam(bam_file, threads):
    command = f"samtools index {bam_file} -@ {threads}"
    print(f"Indexing BAM file with command: {command}")
    subprocess.run(command, shell=True, check=True)

# Function to merge multiple BAM files
def merge_bams(bam_files, output_bam, threads):
    command = f"samtools merge -@ {threads} {output_bam} " + " ".join(bam_files)
    print(f"Merging BAM files with command: {command}")
    subprocess.run(command, shell=True, check=True)
    return output_bam

# Function to get multiple files
def find_fastq_files(output_dir, prefix):
    patterns = [
        f"{prefix}_*.fastq",
        f"{prefix}_*.fastq.gz",
        f"{prefix}_*.fq",
        f"{prefix}_*.fq.gz",
    ]

    files = []
    for pattern in patterns:
        files.extend(glob.glob(os.path.join(output_dir, pattern)))

    return sorted(files)

# Remove FutureWarnings
warnings.simplefilter(action='ignore', category=FutureWarning)

# Main function
def main(reference_genome, modified_genome, method_file, method, coverage, allele_frequency, output_dir, technology, threads):
    print(f'Reference genome: {reference_genome}')
    print(f'Modified genome: {modified_genome}')
    print(f'Method file: {method_file}')
    print(f'Method: {method}')
    print(f'Coverage: {coverage}')
    print(f'Allele frequency: {allele_frequency}')
    print(f'Technology: {technology}')
    print(f'Threads: {threads}')
    print(f'Output directory: {output_dir}')

    # Calculate coverage of modified and reference genome
    reference_coverage = int(coverage * (1 - allele_frequency))
    modified_coverage = int(coverage * allele_frequency)

    # Create the output directory
    os.makedirs(output_dir, exist_ok=True)

    # Generate reads for reference and modified genomes using PBSIM
    # Generate reads depending on allele frequency
    if allele_frequency == 0:
        run_pbsim(reference_genome, method_file, method, coverage, output_dir, 'Reference_reads')

    elif allele_frequency == 1:
        run_pbsim(modified_genome, method_file, method, coverage, output_dir, 'Modified_reads')

    else:
        run_pbsim(reference_genome, method_file, method, reference_coverage, output_dir, 'Reference_reads')
        run_pbsim(modified_genome, method_file, method, modified_coverage, output_dir, 'Modified_reads')

    # Search all FASTQ files for reference and modified genomes
    reference_fastq_files = find_fastq_files(output_dir, "Reference_reads")
    modified_fastq_files = find_fastq_files(output_dir, "Modified_reads")

    bam_files = []

    # CASE 1: Only reference
    if allele_frequency == 0:
        for i, ref_fastq in enumerate(reference_fastq_files):
            output_bam = os.path.join(output_dir, f"combined_alignment_{i + 1}.bam")
            run_minimap2(reference_genome, ref_fastq, "", output_bam, technology, threads)

            sorted_bam_file = sort_bam(output_bam, threads)
            index_bam(sorted_bam_file, threads)
            bam_files.append(sorted_bam_file)

    # CASE 2: Only modified
    elif allele_frequency == 1:
        for i, mod_fastq in enumerate(modified_fastq_files):
            output_bam = os.path.join(output_dir, f"combined_alignment_{i + 1}.bam")
            run_minimap2(reference_genome, mod_fastq, "", output_bam, technology, threads)

            sorted_bam_file = sort_bam(output_bam, threads)
            index_bam(sorted_bam_file, threads)
            bam_files.append(sorted_bam_file)

    # CASE 3: Both
    else:
        if len(reference_fastq_files) != len(modified_fastq_files):
            raise ValueError("The number of reference and modified FASTQ files do not match.")

        for i in range(len(reference_fastq_files)):
            output_bam = os.path.join(output_dir, f"combined_alignment_{i + 1}.bam")

            run_minimap2(
                reference_genome,
                reference_fastq_files[i],
                modified_fastq_files[i],
                output_bam,
                technology,
                threads
            )

            sorted_bam_file = sort_bam(output_bam, threads)
            index_bam(sorted_bam_file, threads)
            bam_files.append(sorted_bam_file)

    # Merge all BAM files into a single BAM file
    final_bam_file = os.path.join(output_dir, 'combined_final_alignment.bam')
    final_bam = merge_bams(bam_files, final_bam_file, threads)

    # Sort and index the final merged BAM file
    sorted_bam_file = sort_bam(final_bam, threads)
    index_bam(sorted_bam_file, threads)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Generate and align reads from reference and modified genome')
    parser.add_argument('--reference_genome', type=str, required=False, help='Path to the reference genome (.FASTA)')
    parser.add_argument('--modified_genome', type=str, required=False, help='Path to the modified genome (.FASTA)')
    parser.add_argument('--method_file', type=str, required=True, help='Path to the method file (.MODEL or .FASTQ)')
    parser.add_argument('--method', type=str, required=True, choices=['quality_score', 'error_model', 'training'], help='Method to use with PBSIM (quality_score,error_model,training)')
    parser.add_argument('--coverage', type=int, required=True, help='Coverage (100 for 100x)')
    parser.add_argument('--allele_frequency', type=float, required=True, help='Allele frequency (0.25 for 25%)')
    parser.add_argument('--output_dir', type=str, required=True, help='Name of the parent output directory where results will be saved')
    parser.add_argument('--technology', type=str, required=True, choices=['ONT', 'PB', 'HiFi'], help='Sequencing technology to use (ONT, PB, HiFi)')
    parser.add_argument('--threads', type=int, default=1, required=False, help='Number of threads to use for Minimap2, and Samtools')

    args = parser.parse_args()

    # Validate genome inputs depending on allele frequency
    if args.allele_frequency == 0:
        if not args.reference_genome:
            parser.error("--reference_genome is required when allele_frequency = 0")

    elif args.allele_frequency == 1:
        if not args.modified_genome:
            parser.error("--modified_genome is required when allele_frequency = 1")

    else:
        if not args.reference_genome or not args.modified_genome:
            parser.error("--reference_genome and --modified_genome are required when 0 < allele_frequency < 1")

    main(args.reference_genome, args.modified_genome, args.method_file, args.method, args.coverage, args.allele_frequency, args.output_dir, args.technology, args.threads)
    
    method_file = os.path.abspath(method_file)
    output_prefix = os.path.join(output_dir, output_reference)

    if method == 'quality_score':
        command = f"pbsim --strategy wgs --method qshmm --qshmm {method_file} --depth {depth} --genome {genome} --prefix {output_prefix}"
    elif method == 'error_model':
        command = f"pbsim --strategy wgs --method errhmm --errhmm {method_file} --depth {depth} --genome {genome} --prefix {output_prefix}"
    elif method == 'training':
        command = f"pbsim --strategy wgs --method sample --sample {method_file} --depth {depth} --genome {genome} --prefix {output_prefix}"
    else:
        raise ValueError(f"Unknown method: {method}")

    print(f"Running PBSIM with command: {command}")
    subprocess.run(command, shell=True, check=True)

    return output_prefix


# Function to align reads using Minimap2
def run_minimap2(reference_file, fastq_file_1, fastq_file_2, output_bam, technology, threads):
    if technology == 'ONT':
        command = f"minimap2 -ax map-ont {reference_file} {fastq_file_1} {fastq_file_2} -t {threads}"
    elif technology == 'PB':
        command = f"minimap2 -ax map-pb {reference_file} {fastq_file_1} {fastq_file_2} -t {threads}"
    elif technology == 'HiFi':
        command = f"minimap2 -ax map-hifi {reference_file} {fastq_file_1} {fastq_file_2} -t {threads}"
    else:
        raise ValueError(f"Unknown technology: {technology}")

    command += f" | samtools view -bS -o {output_bam} -@ {threads}"
    print(f"Running Minimap2 with command: {command}")
    subprocess.run(command, shell=True, check=True)

# Function to sort BAM file
def sort_bam(bam_file, threads):
    sorted_bam_file = bam_file.replace('.bam', '.sorted.bam')
    command = f"samtools sort {bam_file} -o {sorted_bam_file} -@ {threads}"
    print(f"Sorting BAM file with command: {command}")
    subprocess.run(command, shell=True, check=True)
    return sorted_bam_file

# Function to index BAM file
def index_bam(bam_file, threads):
    command = f"samtools index {bam_file} -@ {threads}"
    print(f"Indexing BAM file with command: {command}")
    subprocess.run(command, shell=True, check=True)

# Function to merge multiple BAM files
def merge_bams(bam_files, output_bam, threads):
    command = f"samtools merge -@ {threads} {output_bam} " + " ".join(bam_files)
    print(f"Merging BAM files with command: {command}")
    subprocess.run(command, shell=True, check=True)
    return output_bam

# Remove FutureWarnings
warnings.simplefilter(action='ignore', category=FutureWarning)

# Main function
def main(reference_genome, modified_genome, method_file, method, coverage, allele_frequency, output_dir, technology, threads):
    print(f'Reference genome: {reference_genome}')
    print(f'Modified genome: {modified_genome}')
    print(f'Method file: {method_file}')
    print(f'Method: {method}')
    print(f'Coverage: {coverage}')
    print(f'Allele frequency: {allele_frequency}')
    print(f'Technology: {technology}')
    print(f'Threads: {threads}')
    print(f'Output directory: {output_dir}')

    # Calculate coverage of modified and reference genome
    reference_coverage = int(coverage * (1 - allele_frequency))
    modified_coverage = int(coverage * allele_frequency)

    # Create the output directory
    os.makedirs(output_dir, exist_ok=True)

    # Generate reads for reference and modified genomes using PBSIM
    run_pbsim(reference_genome, method_file, method, reference_coverage, output_dir, 'Reference_reads')
    run_pbsim(modified_genome, method_file, method, modified_coverage, output_dir, 'Modified_reads')
    
    # Search all FASTQ files for reference and modified genomes
    reference_fastq_files = sorted(glob.glob(os.path.join(output_dir, "Reference_reads_*.fastq")))
    modified_fastq_files = sorted(glob.glob(os.path.join(output_dir, "Modified_reads_*.fastq")))

    if len(reference_fastq_files) != len(modified_fastq_files):
        raise ValueError("The number of reference and modified FASTQ files do not match.")

    # Empty list for BAM files
    bam_files = []

    # Minimap + Sort + Index for each sample
    for i in range(len(reference_fastq_files)):
        reference_fastq_file = reference_fastq_files[i]
        modified_fastq_file = modified_fastq_files[i]
        output_bam = os.path.join(output_dir, f"combined_alignment_{i + 1}.bam")

        # Run Minimap2 for each pair of files
        run_minimap2(reference_genome, reference_fastq_file, modified_fastq_file, output_bam, technology, threads)
        
        # Sort and index the generated BAM file
        sorted_bam_file = sort_bam(output_bam, threads)
        index_bam(sorted_bam_file, threads)

        # Add the sorted BAM to the list of BAM files to merge later
        bam_files.append(sorted_bam_file)

    # Merge all BAM files into a single BAM file
    final_bam_file = os.path.join(output_dir, 'combined_final_alignment.bam')
    final_bam = merge_bams(bam_files, final_bam_file, threads)

    # Sort and index the final merged BAM file
    sorted_bam_file = sort_bam(final_bam, threads)
    index_bam(sorted_bam_file, threads)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Generate and align reads from reference and modified genome')
    parser.add_argument('--reference_genome', type=str, required=True, help='Path to the reference genome (.FASTA)')
    parser.add_argument('--modified_genome', type=str, required=True, help='Path to the modified genome (.FASTA)')
    parser.add_argument('--method_file', type=str, required=True, help='Path to the method file (.MODEL or .FASTQ)')
    parser.add_argument('--method', type=str, required=True, choices=['quality_score', 'error_model', 'training'], help='Method to use with PBSIM (quality_score,error_model,training)')
    parser.add_argument('--coverage', type=int, required=True, help='Coverage (100 for 100x)')
    parser.add_argument('--allele_frequency', type=float, required=True, help='Allele frequency (0.25 for 25%)')
    parser.add_argument('--output_dir', type=str, required=True, help='Name of the parent output directory where results will be saved')
    parser.add_argument('--technology', type=str, required=True, choices=['ONT', 'PB', 'HiFi'], help='Sequencing technology to use (ONT, PB, HiFi)')
    parser.add_argument('--threads', type=int, default=1, required=False, help='Number of threads to use for Minimap2, and Samtools')

    args = parser.parse_args()
    main(args.reference_genome, args.modified_genome, args.method_file, args.method, args.coverage, args.allele_frequency, args.output_dir, args.technology, args.threads)
