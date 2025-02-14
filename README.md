# SVModeller

SVModeller is a computational tool to simulate synthetic human haplotypes containing embedded structural variants (SV). This simulator has been trained using an extensive catalogue of curated polymorphic SVs identified in a dataset comprising 1.019 samples from the 1000 Genomes Project, sequenced with Oxford Nanopore long-read technology (Schloissnig et al., 2024). This dataset provides detailed sequence information and annotation of SV classes that has not been previously available. 

One of the key innovations of this simulator is its ability to include large repetitive insertions, such as Variable Number of Tandem Repeats (VNTR), mobile element insertions (MEI), transduction-related insertions, and mitochondrial insertions; all of them taking into consideration their specific features such as polyadenylation tails, SVA hexamers or TSDs, between others, addressing the limitations of prior simulators. This computational tool will enable the research community to create gold-standard sets of synthetic human genomes, which will be crucial for benchmarking and evaluating methods for detecting and annotating SVs using long-read sequencing data.

The simulator is implemented in python (Van Rossum & Drake, 2009), and is structured in 5 different independent modules: 
- Module 1: Process insertion features distributions from a variant calling file (VCF).
- Module 2: Generate insertion events from user-defined distributions or those derived from Module 1.
- Module 3: Generate deletion events sampling data stored in a VCF.
- Module 4: Modify a reference genome.
- Module 5: Generate sub-clonal variant reads at different coverage and allele frequency levels.

Structural Variants generated: Mobile Element Insertions (MEIs) of Alu, LINE-1 and SVA, MEIs mediating transductions, mitochondrial insertions (NUMT), VNTRs, duplications, inversions, deletions, inverted duplications, orphan insertions, and deletions.

## Modules description, input and output
<img width="881" alt="Figure_SVModeller" src="https://github.com/user-attachments/assets/3812ce77-67a2-4721-820a-fadec290bbb6" />

### Module 1
This module aims to process the insertion data from the VCF obtaining different distributions such as how the different events are genome-wide spread, the distribution of the different features that constitute each event, and the probability of finding each event on the genome.

**Input:**
- VCF with insertion data
- Chromosomes length (.txt)
- Window size for genome segmentation (by default 1 Mega base)
  
**Output:**
- Genome-wide distribution (.tsv)
- Insertion features (.tsv)
- Event probabilities (.tsv)

### Module 2
Second module objective is to generate from determine distributions the final sequences that will be inserted in a genome. It takes the genome-wide and insertion features distributions, even from module 1 or user-defined ones with the same structure as the ones derived from module 1. From these distributions samples different values to build the new events and returns their generated sequence.

**Input:**
- Genome-Wide Distribution (.tsv)
- Insertion Features (.tsv)
- Event Probabilities (.tsv)
- Total number of events to simulate 
- Table with source loci to LINE-1 transductions (.tsv)
- Table with source loci to SVA transductions (.tsv)
- Chromosomes length (.txt)
- Consensus sequences (.fasta)
- Reference genome (.fasta)

**Output:**
- New insertion events sequences with their corresponding features (.tsv)

### Module 3
Third module generates deletions events. It takes the data from VCF and based on the determine number of events to simulate, selects the regions to be deleted, and returns a data frame with the regions to be removed from the genome.

**Input:**
- VCF SVs classified deletions
- VCF SVs unclassified deletions
- Total number of events to simulate
- Chromosomes length (.txt)
- Window size for genome segmentation (by default 1 Mega base)

**Output:**
- Deletion regions (.tsv)

### Module 4
Fourth module results in a modified genome. It takes data frames of insertion and deletion data (even from modules 2 and 3 or user-defined) and modifies a reference genome with the determine information. The result is the reference genome with the specified changes.

**Input:**
- Insertion events data (.tsv)
- Deletion events data (.tsv)
- Reference genome (.fasta)

**Output:**
- Modified reference genome (.fasta)

### Module 5
Fifth module aims to generate sub-clonal variant reads with different coverage and allele frequency levels. It takes the reference and modified genomes, as well as the method file (based on error model, quality score, or fastQ files) and the corresponding method file, the technology of the simulated reads (Oxford Nanopore - ONT, PacBio HiFi, PacBio HiFi), and desired coverage and allele frequency.

**Input:**
- Reference genome (.fasta)
- Modified reference genome (.fasta)
- Method to generate reads
- Method file
- Allele frequency
- Coverage
- Output directory
- Technology of generated reads

**Output:**
- Alignment (.bam)
- Modified genome reads (.fastq)
- Reference genome reads (.fastq)


## Data
Available data to run the modules available at: https://zenodo.org/records/14629433?token=eyJhbGciOiJIUzUxMiJ9.eyJpZCI6IjkxODlmMGVmLTU4ZmMtNGQ5NC04MmRmLWNiNTc3NWRlZGEzYSIsImRhdGEiOnt9LCJyYW5kb20iOiIzMTQyNDg4MjcyYTI2ZDRmZTI3MTcwOGJkNTEzY2RiNyJ9.0ikZdwVU-K6ffwjmEb3HudqHvyz52Umt0XZyO91t_1ffEeP3VTV9iKaRhTWBf8ABoF03R24WMKzMu23yDjpX3Q 

## Copyright
Copyright © 2025 Ismael Vera Muñoz. All rights reserved.
