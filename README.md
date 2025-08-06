## make-the-snake
A Snakemake-based pipeline for automated variant calling and annotation from raw sequencing reads using publicly available SRA and GenBank data.

### Overview
This pipeline processes raw FASTQ reads through quality control, alignment, variant calling, filtering, and annotation. It is designed for small-scale data and instructional use, with a focus on reproducibility and modularity.

### Tools Used
- Snakemake
- BWA
- samtools
- GATK
- SnpEff
- FastQC
- NCBI E-utilities (efetch, prefetch, fastq-dump)

### Input
SRA accession number (default: SRR1972739)
NCBI reference genome ID (default: AF086833.2)

### Workflow Steps
Download reference genome and sequencing reads
Perform quality control on FASTQ files
Index the reference genome (samtools, BWA, GATK)
Align reads to reference using BWA
Convert SAM to sorted and deduplicated BAM
Call variants with GATK HaplotypeCaller
Filter variants based on quality metrics
Download GenBank file and build custom SnpEff database
Annotate filtered variants with SnpEff

### Output
Sorted and indexed BAM files
Raw and filtered VCF files
Annotated VCF file
SnpEff HTML summary report
QC reports from FastQC
