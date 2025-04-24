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
- VCF with insertion data _(VCF_Insertions.vcf)_
- Chromosomes length _(chr_length.txt)_
- OPTIONAL: Window size for genome segmentation, by default 1 Mega base _(integer number)_
  
**Output:**
- Genome-wide distribution _(Genome_Wide_Distribution.tsv.tsv)_
- Insertion features _(Insertion_Features.tsv)_
- Event probabilities _(Probabilities.tsv)_
- List of VNTR motifs _(Separated_Motifs.tsv)_
- List of SVA VNTR motifs _(SVA_VNTR_Motifs.txt)_

### Module 2
Second module objective is to generate from determine distributions the final sequences that will be inserted in a genome. It takes the genome-wide and insertion features distributions, even from module 1 or user-defined ones with the same structure as the ones derived from module 1. From these distributions samples different values to build the new events and returns their generated sequence.

**Input:**
- Genome-Wide Distribution _(Genome_Wide_Distribution.tsv)_
- Insertion Features _(Insertion_Features.tsv)_
- Event Probabilities or Number of each event to simulate _(Probabilities.tsv)_
- OPTIONAL: Just in case of providing probabilities, total number of events to simulate _(integer number)_
- Table with source loci to LINE-1 transductions _(source_loci_LINE1.tsv)_
- Table with source loci to SVA transductions _(source_loci_SVA.tsv)_
- Consensus sequences _(consensus_sequences_complete.fa)_
- Reference genome _(chm13v2.0.fa)_
- List of VNTR motifs _(Separated_Motifs.tsv)_
- List of SVA VNTR motifs _(SVA_VNTR_Motifs.txt)_

**Output:**
- New insertion events sequences with their corresponding features _(Insertions_table.tsv)_
- OPTIONAL: Variant Calling File (VCF) with insertion data

**Mandatory input columns**
- Genome-Wide Distribution _(Genome_Wide_Distribution.tsv.tsv)_: 'window', 'beg', 'end', each event conformation as a column (i.e. 'Alu__FOR+POLYA', 'VNTR', 'INV_DUP')
- Insertion features _(Insertion_Features.tsv)_: 'Event', 'Length', 'Strand', 'TSD_Length', 'TD_5', 'TD_3', 'SVA_Hexamer', 'SVA_VNTR_Length', 'TD_orphan_Length', 'VNTR_Num_Motifs', 'PolyA_Length_1', 'PolyA_Length_2', 'FOR', 'TRUN', 'REV', 'DEL', 'DUP'
- Event Probabilities _(Probabilities.tsv)_: 'Event' and 'Probability' (probabilities between 0 and 1)or 'Number' (number of each event to simulate)
   
### Module 3
Third module generates deletions events. It takes the data from VCF and based on the determine number of events to simulate, selects the regions to be deleted, and returns a data frame with the regions to be removed from the genome.

**Input:**
- VCF with deletion data _(VCF_Deletions.vcf)_
- Total number of events to simulate _(integer number)_
- Chromosomes length _(chr_length.txt)_
- OPTIONAL: Reference genome _(chm13v2.0.fa)_ just if VCF file is desired
- OPTIONAL: Window size for genome segmentation, by default 1 Mega base _(integer number)_

**Output:**
- Deletion regions _(Deletions_table.tsv)_
- OPTIONAL: Variant Calling File (VCF) with deletion data

### Module 4
Fourth module results in a modified genome. It takes data frames of insertion and deletion data (even from modules 2 and 3 or user-defined) and modifies a reference genome with the determine information. The result is the reference genome with the specified changes.

**Input:**
- Events to mofidy reference genome (.tsv)
- Reference genome (.fasta)
- OPTIONAL: additional events table (.tsv)

**Output:**
- Modified reference genome _(Modified_Reference_Genome.fasta)_
- Table with added events and their current positions in the genome _(Sorted_Genomic_Events.tsv)_

**Mandatory input columns**
- Events to mofidy reference genome: '#ref', 'beg', 'Length', 'Event_Type', 'Sequence_Insertion'

### Module 5
Fifth module aims to generate sub-clonal variant reads with different coverage and allele frequency levels. It takes the reference and modified genomes, as well as the method file (based on error model, quality score, or fastQ files) and the corresponding method file, the technology of the simulated reads (Oxford Nanopore - ONT, PacBio - PB, PacBio-HiFi - HiFi), and desired coverage and allele frequency. Please ensure consistency in the technology used across all files

**Input:**
- Reference genome (.fasta)
- Modified reference genome (.fasta)
- Method to generate reads _(choices='quality_score', 'error_model', 'training')_
- Method file _(string, '.model' or '.fastq' i.e. ERRHMM-ONT-HQ.model)_
- Allele frequency _(decimal number between 0-1. i.e. 0.25)_
- Coverage _(integer number: 30 for 30x)_
- Output directory _(file name)_
- Technology of generated reads _(ONT, PB, HiFi)_

**Output:**
- Alignment _(combined_final_alignment.bam)_
- Modified genome reads (.fastq)
- Reference genome reads (.fastq)


## Data
Available data to run the modules for versions 0.1 and 0.2 available at:

- Version 0.2:
[Link to Zenodo](https://zenodo.org/records/14967709?token=eyJhbGciOiJIUzUxMiJ9.eyJpZCI6IjQ0NDU5OWVjLTZiM2EtNDllOC05YzNkLTE3YjNiM2NhNjVmNCIsImRhdGEiOnt9LCJyYW5kb20iOiI2NGFkOGIxNGY4YmQyOWFlNzlhZGU1NzEwNzZiYmEzNSJ9.ueQ38ZbCp-OBSGps4-_Q5KUDtJ2JoGPhUApvrzHgtZ6j5IqzLIdrik4BEPhVqteQl4yxZ9en_A57mmS_U3zayQ)

- Version 0.1:
[Link to Zenodo](https://zenodo.org/records/14629433?token=eyJhbGciOiJIUzUxMiJ9.eyJpZCI6IjkxODlmMGVmLTU4ZmMtNGQ5NC04MmRmLWNiNTc3NWRlZGEzYSIsImRhdGEiOnt9LCJyYW5kb20iOiIzMTQyNDg4MjcyYTI2ZDRmZTI3MTcwOGJkNTEzY2RiNyJ9.0ikZdwVU-K6ffwjmEb3HudqHvyz52Umt0XZyO91t_1ffEeP3VTV9iKaRhTWBf8ABoF03R24WMKzMu23yDjpX3Q)


## Additional scripts
- Filter_VCF_information.py: to remove from VCF files the information field, keeping just the event length.


## Developers
SVModeller has been developed by Ismael Vera Muñoz at the Repetitive DNA Biology (REPBIO) Lab at the Centre for Genomic Regulation (CRG) (2024-2025).

## License
SVModeller is distributed under the AGPL-3.0. Consult the [LICENSE](https://github.com/REPBIO-LAB/SVModeller/blob/main/LICENSE) file for more information.
