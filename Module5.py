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
import random
import numpy as np
import tempfile

def set_seed(seed: int = 42):
    """Sets random seeds across all libraries used in SVModeller

    to ensure reproducible results.
    """
    if seed is None:
        return

    # Set standard Python random seed
    random.seed(seed)

    # Set NumPy random seed
    np.random.seed(seed)

    # Set Python hash seed for consistent dictionary/set ordering
    os.environ['PYTHONHASHSEED'] = str(seed)

# Function to run PBSIM to generate synthetic reads for reference and modified genomes
def run_pbsim(genome, method_file, method, depth, output_dir, output_reference, seed):
    if depth == 0:
        print("Depth is 0. Skipping PBSIM execution.")
        return
    
    method_file = os.path.abspath(method_file)
    output_prefix = os.path.join(output_dir, output_reference)

    if method == 'quality_score':
        command = f"pbsim --strategy wgs --method qshmm --qshmm {method_file} --depth {depth} --genome {genome} --prefix {output_prefix} --seed {seed}"
    elif method == 'error_model':
        command = f"pbsim --strategy wgs --method errhmm --errhmm {method_file} --depth {depth} --genome {genome} --prefix {output_prefix} --seed {seed}"
    elif method == 'training':
        command = f"pbsim --strategy wgs --method sample --sample {method_file} --depth {depth} --genome {genome} --prefix {output_prefix} --seed {seed}"
    else:
        raise ValueError(f"Unknown method: {method}")

    print(f"Running PBSIM with command: {command}")
    subprocess.run(command, shell=True, check=True)

    return output_prefix

def run_insilicoseq_single_genome(genome, short_read_model, depth, output_dir, output_reference, seed):
    if depth == 0:
        print("Depth is 0. Skipping InsilicoSeq execution.")
        return

    if short_read_model not in ["HiSeq", "MiSeq", "MiSeq-20", "MiSeq-24", "MiSeq-28", "MiSeq-32", "MiSeq-36", "NextSeq", "NovaSeq"]:
        raise ValueError(f"Unknown simulation model for InsilicoSeq: {short_read_model}")

    output_prefix = os.path.join(output_dir, output_reference)

    coverage_content = f"{genome}\t{depth}\n"
    
    with tempfile.NamedTemporaryFile(mode="w+", delete=True) as cov_file:
        cov_file.write(coverage_content)
        cov_file.flush()

        cmd = [
            "iss",
            "generate",
            "--draft",
            str(genome),
            "--coverage_file",
            cov_file.name,
            "--model",
            str(short_read_model),
            "--seed",
            str(seed),
            "--output",
            str(output_prefix),
        ]

        print(f"Running InsilicoSeq with command: {' '.join(cmd)}")
        subprocess.run(cmd, check=True)

    return output_prefix

def run_insilicoseq_diploid(hap1, hap2, short_read_model, depth, output_dir, output_reference, seed):
    if depth == 0:
        print("Depth is 0. Skipping InsilicoSeq execution.")
        return

    if short_read_model not in ["HiSeq", "MiSeq", "MiSeq-20", "MiSeq-24", "MiSeq-28", "MiSeq-32", "MiSeq-36", "NextSeq", "NovaSeq"]:
        raise ValueError(f"Unknown simulation model for InsilicoSeq: {short_read_model}")

    output_prefix = os.path.join(output_dir, output_reference)
    
    diploid_depth = depth/2
    coverage_content = f"{hap1}\t{diploid_depth}\n{hap2}\t{diploid_depth}\n"

    with tempfile.NamedTemporaryFile(mode="w+", delete=True) as cov_file:
        cov_file.write(coverage_content)
        cov_file.flush()

        cmd = [
            "iss",
            "generate",
            "--draft",
            str(hap1),
            str(hap2),
            "--coverage_file",
            cov_file.name,
            "--model",
            str(short_read_model),
            "--seed",
            str(seed),
            "--output",
            str(output_prefix),
        ]

        print(f"Running InsilicoSeq with command: {' '.join(cmd)}")
        subprocess.run(cmd, check=True)

    return output_prefix

# Function to align reads using Minimap2
def run_minimap2(reference_file, fastq_file_1, fastq_file_2, output_bam, technology, threads, seed):
    fastqs = fastq_file_1
    if fastq_file_2:
        fastqs += f" {fastq_file_2}"

    if technology == 'ONT':
        command = f"minimap2 -ax map-ont {reference_file} {fastqs} -t {threads} --seed {seed}"
    elif technology == 'PB':
        command = f"minimap2 -ax map-pb {reference_file} {fastqs} -t {threads} --seed {seed}"
    elif technology == 'HiFi':
        command = f"minimap2 -ax map-hifi {reference_file} {fastqs} -t {threads} --seed {seed}"
    else:
        raise ValueError(f"Unknown technology: {technology}")

    command += f" | samtools view -bS -o {output_bam} -@ {threads}"
    subprocess.run(command, shell=True, check=True)

def run_bwa_short_reads(reference_file, fastq_file_1, fastq_file_2, output_bam, threads):
    #the random seed is not needed in bwa, however the multi-thread change might introduce a bit difference if using different threads for re-produce
    command = f"bwa mem -t {threads} -R '@RG\\tID:{output_bam}\\tSM:{output_bam}\\tPL:ILLUMINA' {reference_file} {fastq_file_1} {fastq_file_2}"

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

def find_paired_fastq_files(output_dir, prefix):
    patterns_1 = [
        f"{prefix}_R1.fastq",
        f"{prefix}_R1.fastq.gz",
        f"{prefix}_R1.fq",
        f"{prefix}_R1.fq.gz",
    ]

    patterns_2 = [
        f"{prefix}_R2.fastq",
        f"{prefix}_R2.fastq.gz",
        f"{prefix}_R2.fq",
        f"{prefix}_R2.fq.gz",
    ]

    files_1 = []
    files_2 = []
    for pattern in patterns_1:
        files_1.extend(glob.glob(os.path.join(output_dir, pattern)))
    for pattern in patterns_2:
        files_2.extend(glob.glob(os.path.join(output_dir, pattern)))

    return sorted(files_1), sorted(files_2)

# Remove FutureWarnings
warnings.simplefilter(action='ignore', category=FutureWarning)

# Main function
def main(reference_genome, modified_genome, method_file, method, coverage, allele_frequency, output_dir, technology, short_read_model, mode, threads, seed):
    # argparse's append action returns a list; single-genome simulators need a path.
    if mode != 'diploid':
        if isinstance(modified_genome, list):
            if len(modified_genome) != 1:
                raise ValueError("Single-genome mode requires exactly one modified genome path.")
            modified_genome = modified_genome[0]

    print(f'Reference genome: {reference_genome}')
    print(f'Modified genome: {modified_genome}')
    print(f'Method file: {method_file}')
    print(f'Method: {method}')
    print(f'Coverage: {coverage}')
    print(f'Allele frequency: {allele_frequency}')
    print(f'Technology: {technology}')
    print(f'Threads: {threads}')
    print(f'Output directory: {output_dir}')

    # Set seed
    set_seed(seed)

    # Calculate coverage of modified and reference genome
    if allele_frequency is not None:
        reference_coverage = int(coverage * (1 - allele_frequency))
        modified_coverage = int(coverage * allele_frequency)

    # Create the output directory
    os.makedirs(output_dir, exist_ok=True)

    # Generate reads for reference and modified genomes using PBSIM
    # Generate reads depending on allele frequency
    if technology == 'Illumina':
        if mode == 'haploid':
            if allele_frequency == 0:
                run_insilicoseq_single_genome(reference_genome, short_read_model, coverage, output_dir, 'Reference_reads', seed)
            elif allele_frequency == 1:
                run_insilicoseq_single_genome(modified_genome, short_read_model, coverage, output_dir, 'Modified_reads', seed)
        
            else:
                run_insilicoseq_single_genome(reference_genome, short_read_model, reference_coverage, output_dir, 'Reference_reads', seed)
                run_insilicoseq_single_genome(modified_genome, short_read_model, modified_coverage, output_dir, 'Modified_reads', seed + 1)
        elif mode == 'diploid':
            hap1 = modified_genome[0]
            hap2 = modified_genome[1]
            run_insilicoseq_diploid(hap1, hap2, short_read_model, coverage, output_dir, 'Modified_reads', seed)
    else:
        if mode == 'haploid':
            if allele_frequency == 0:
                run_pbsim(reference_genome, method_file, method, coverage, output_dir, 'Reference_reads', seed)
    
            elif allele_frequency == 1:
                run_pbsim(modified_genome, method_file, method, coverage, output_dir, 'Modified_reads', seed)
    
            else:
                run_pbsim(reference_genome, method_file, method, reference_coverage, output_dir, 'Reference_reads', seed)
                run_pbsim(modified_genome, method_file, method, modified_coverage, output_dir, 'Modified_reads', seed + 1)
        elif mode == 'diploid':
            hap1 = modified_genome[0]
            hap2 = modified_genome[1]
            run_pbsim(hap1, method_file, method, coverage/2, output_dir, 'Modified_reads_1', seed)
            run_pbsim(hap2, method_file, method, coverage/2, output_dir, 'Modified_reads_2', seed + 1)

    bam_files = []

    if technology == "Illumina":
        # Search all FASTQ files for reference and modified genomes
        reference_fastq_files_1, reference_fastq_files_2= find_paired_fastq_files(output_dir, "Reference_reads")
        modified_fastq_files_1, modified_fastq_files_2= find_paired_fastq_files(output_dir, "Modified_reads")

        if len(reference_fastq_files_1) != len(reference_fastq_files_2):
            raise ValueError("The number of reference R1 and R2 FASTQ files do not match.")

        if len(modified_fastq_files_1) != len(modified_fastq_files_2):
            raise ValueError("The number of modified R1 and R2 FASTQ files do not match.")
        
        if mode == "diploid":
            for i, (fq1, fq2) in enumerate(zip(modified_fastq_files_1, modified_fastq_files_2)):
                output_bam = os.path.join(output_dir, f"combined_alignment_{i + 1}.bam")
                run_bwa_short_reads(reference_genome, fq1, fq2, output_bam, threads)
            
                sorted_bam_file = sort_bam(output_bam, threads)
                index_bam(sorted_bam_file, threads)
                bam_files.append(sorted_bam_file)

        elif mode == "haploid":
            # CASE 1: Only reference
            if allele_frequency == 0:
                for i, (fq1, fq2) in enumerate(zip(reference_fastq_files_1, reference_fastq_files_2)):
                    output_bam = os.path.join(output_dir, f"combined_alignment_{i + 1}.bam")
                    run_bwa_short_reads(reference_genome, fq1, fq2, output_bam, threads)
    
                    sorted_bam_file = sort_bam(output_bam, threads)
                    index_bam(sorted_bam_file, threads)
                    bam_files.append(sorted_bam_file)
    
            # CASE 2: Only modified
            elif allele_frequency == 1:
                for i, (fq1, fq2) in enumerate(zip(modified_fastq_files_1, modified_fastq_files_2)):
                    output_bam = os.path.join(output_dir, f"combined_alignment_{i + 1}.bam")
                    run_bwa_short_reads(reference_genome, fq1, fq2, output_bam, threads)
    
                    sorted_bam_file = sort_bam(output_bam, threads)
                    index_bam(sorted_bam_file, threads)
                    bam_files.append(sorted_bam_file)
    
            # CASE 3: Both
            else:
                if len(reference_fastq_files_1) != len(modified_fastq_files_1):
                    raise ValueError("The number of reference and modified FASTQ files do not match.")
    
                for i in range(len(reference_fastq_files_1)):
                    output_bam = os.path.join(output_dir, f"reference_alignment_{i + 1}.bam")
                    run_bwa_short_reads(reference_genome, reference_fastq_files_1[i], reference_fastq_files_2[i], output_bam, threads)
    
                    sorted_bam_file = sort_bam(output_bam, threads)
                    index_bam(sorted_bam_file, threads)

                    bam_files.append(sorted_bam_file)

                    output_bam = os.path.join(output_dir, f"modified_alignment_{i + 1}.bam")
                    run_bwa_short_reads(reference_genome, modified_fastq_files_1[i], modified_fastq_files_2[i], output_bam, threads)
    
                    sorted_bam_file = sort_bam(output_bam, threads)
                    index_bam(sorted_bam_file, threads)

                    bam_files.append(sorted_bam_file)
        
    else:
        if mode == "diploid":
            # Search all FASTQ files for reference and modified genomes
            modified_fastq_files_1 = find_fastq_files(output_dir, "Modified_reads_1")
            modified_fastq_files_2 = find_fastq_files(output_dir, "Modified_reads_2")

            if not modified_fastq_files_1 or not modified_fastq_files_2:
                raise RuntimeError("Missing FASTQ output for one or both haplotypes.")

            modified_fastq_files = modified_fastq_files_1 + modified_fastq_files_2

            for i, mod_fastq in enumerate(modified_fastq_files):
                output_bam = os.path.join(output_dir, f"combined_alignment_{i + 1}.bam")
                run_minimap2(reference_genome, mod_fastq, "", output_bam, technology, threads, seed)
                sorted_bam_file = sort_bam(output_bam, threads)
                index_bam(sorted_bam_file, threads)
                bam_files.append(sorted_bam_file)

        elif mode == "haploid":
            # Search all FASTQ files for reference and modified genomes
            reference_fastq_files = find_fastq_files(output_dir, "Reference_reads")
            modified_fastq_files = find_fastq_files(output_dir, "Modified_reads")
            # CASE 1: Only reference
            if allele_frequency == 0:
                for i, ref_fastq in enumerate(reference_fastq_files):
                    output_bam = os.path.join(output_dir, f"combined_alignment_{i + 1}.bam")
                    run_minimap2(reference_genome, ref_fastq, "", output_bam, technology, threads, seed)
    
                    sorted_bam_file = sort_bam(output_bam, threads)
                    index_bam(sorted_bam_file, threads)
                    bam_files.append(sorted_bam_file)
    
            # CASE 2: Only modified
            elif allele_frequency == 1:
                for i, mod_fastq in enumerate(modified_fastq_files):
                    output_bam = os.path.join(output_dir, f"combined_alignment_{i + 1}.bam")
                    run_minimap2(reference_genome, mod_fastq, "", output_bam, technology, threads, seed)
    
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
                        threads,
                        seed
                    )
    
                    sorted_bam_file = sort_bam(output_bam, threads)
                    index_bam(sorted_bam_file, threads)
                    bam_files.append(sorted_bam_file)

    if len(bam_files) < 1:
        raise RuntimeError("No bam file after alignment process")
    # Merge all BAM files into a single BAM file
    final_bam_file = os.path.join(output_dir, 'combined_final_alignment.bam')
    final_bam = merge_bams(bam_files, final_bam_file, threads)

    # Sort and index the final merged BAM file
    sorted_bam_file = sort_bam(final_bam, threads)
    index_bam(sorted_bam_file, threads)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Generate and align reads from reference and modified genome')
    parser.add_argument('--reference_genome', type=str, required=False, help='Path to the reference genome (.FASTA)')
    parser.add_argument('--modified_genome', action="append", type=str, required=False, help='Path to the modified genome (.FASTA)')
    parser.add_argument('--method_file', type=str, help='Path to the method file (.MODEL or .FASTQ)')
    parser.add_argument('--method', type=str, choices=['quality_score', 'error_model', 'training'], help='Method to use with PBSIM (quality_score,error_model,training)')
    parser.add_argument('--coverage', type=int, required=True, help='Coverage (100 for 100x)')
    parser.add_argument('--allele_frequency', type=float, help='Allele frequency (0.25 for 25%%)')
    parser.add_argument('--output_dir', type=str, required=True, help='Name of the parent output directory where results will be saved')
    parser.add_argument('--technology', type=str, required=True, choices=['ONT', 'PB', 'HiFi', 'Illumina'], help='Sequencing technology to use (ONT, PB, HiFi, Illumina)')
    parser.add_argument('--mode', type=str, default='haploid', choices=['haploid', 'diploid'], help='Haploid or diploid for simulation')
    parser.add_argument('--short_read_model', type=str, choices=['HiSeq', 'MiSeq', 'MiSeq-20', 'MiSeq-24', 'MiSeq-28', 'MiSeq-32', 'MiSeq-36', 'NextSeq', 'NovaSeq'], help='(InsilicoSeq) Illumina platform to use (HiSeq, MiSeq, NovaSeq, NextSeq)')
    parser.add_argument('--threads', type=int, default=1, required=False, help='Number of threads to use for Minimap2, and Samtools')
    parser.add_argument('--seed', type=int, required=False, default=42, help='Random seed for reproducibility (default: 42).')

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

    # Separate long and short mode
    if args.technology in ['ONT', 'PB', 'HiFi']:
        if not args.method_file or not args.method:
            parser.error("--method_file and --method are required when choosing long-read platfroms")
    elif args.technology == 'Illumina':
        if not args.short_read_model or not args.mode:
            parser.error("--short_read_model and --mode are required when choosing short-read platfrom")
    else:
        parser.error("The given --technology is not valid (not in ['ONT', 'PB', 'HiFi', 'Illumina'])")

    if args.allele_frequency is not None:
        if args.allele_frequency > 0:
            if args.mode != "diploid":
                if args.modified_genome is None:
                    parser.error("--modified_genome needs 1 input when AF > 0")
                elif len(args.modified_genome) != 1:
                    parser.error("--modified_genome needs 1 input when AF > 0")
            else:
                parser.error("diploid mode don't need AF")
        else:
            if args.mode == "diploid":
                parser.error("diploid mode don't need AF")

    if args.mode == "diploid":
        if args.modified_genome is not None: 
            if len(args.modified_genome) != 2:
                parser.error("--modified_genome needs 2 inputs when choosing diploid mode")
        else:
            parser.error("--modified_genome needs 2 inputs when choosing diploid mode")
    
    main(args.reference_genome, args.modified_genome, args.method_file, args.method, args.coverage, args.allele_frequency, args.output_dir, args.technology, args.short_read_model, args.mode, args.threads, args.seed)
