SRA="SRR1972739"
REF_ID="AF086833.2"
RESULTS_FOLDER="results"
RAW_DIR=f"{RESULTS_FOLDER}/raw"
ALIGNED_DIR=f"{RESULTS_FOLDER}/aligned"
VARIANT_DIR=f"{RESULTS_FOLDER}/variants"
ANNOTATED_DIR=f"{RESULTS_FOLDER}/annotated"
QC_DIR=f"{RESULTS_FOLDER}/qc"
SNPEFF_DIR=f"{RESULTS_FOLDER}/snpEff"
SNPEFF_DATA_DIR=f"{SNPEFF_DIR}/data/reference_db"


rule all:
    input:
        f"{ANNOTATED_DIR}/annotated_variants.vcf",
        f"{SNPEFF_DIR}/snpEff.html"

rule fasta:
    output: f"{RAW_DIR}/reference.fasta"
    shell: 
        """
        mkdir -p {RAW_DIR} {ALIGNED_DIR} {VARIANT_DIR} {ANNOTATED_DIR} {QC_DIR} {SNPEFF_DATA_DIR} {SNPEFF_DIR}
        efetch -db nucleotide -id {REF_ID} -format fasta > {output}
        """

rule data:
    output: f"{RAW_DIR}/{SRA}.fastq"
    shell:
        """
        prefetch {SRA} -O {RAW_DIR}
        fastq-dump -X 10000 {RAW_DIR}/{SRA}/{SRA}.sra -O {RAW_DIR}
        """

rule fastqc:
    input: fastq = f"{RAW_DIR}/{SRA}.fastq"
    output:
        html = f"{QC_DIR}/{SRA}_fastqc.html",
        zip = f"{QC_DIR}/{SRA}_fastqc.zip"
    shell:
        """
        fastqc -o {QC_DIR} {input.fastq}
        """

rule samtools:
    input: fasta = f"{RAW_DIR}/reference.fasta"
    output: f"{RAW_DIR}/reference.fasta.fai"
    shell: "samtools faidx {input.fasta}"

rule bwa:
    input: fasta = f"{RAW_DIR}/reference.fasta"
    output: 
        f"{RAW_DIR}/reference.fasta.amb",
        f"{RAW_DIR}/reference.fasta.ann",
        f"{RAW_DIR}/reference.fasta.bwt",
        f"{RAW_DIR}/reference.fasta.pac",
        f"{RAW_DIR}/reference.fasta.sa"
    shell: "bwa index {input.fasta}"

rule dictionary:
    input: fasta = f"{RAW_DIR}/reference.fasta"
    output: 
        dict = f"{RAW_DIR}/reference.dict",
    shell: "gatk CreateSequenceDictionary -R {input.fasta} -O {output.dict}"

rule alignment:
    input: 
        fasta = f"{RAW_DIR}/reference.fasta",
        fastq = f"{RAW_DIR}/{SRA}.fastq"
    output: 
        sam = f"{ALIGNED_DIR}/aligned.sam"
    shell: 
        "bwa mem -R '@RG\\tID:1\\tLB:lib1\\tPL:illumina\\tPU:unit1\\tSM:sample1' {input.fasta} {input.fastq} > {output.sam}"

rule sam2bam:
    input: sam = f"{ALIGNED_DIR}/aligned.sam"
    output: 
        bam = f"{ALIGNED_DIR}/aligned.sorted.bam"
    shell: 
        "samtools view -b {input.sam} | samtools sort -o {output.bam}"

rule validate_bam:
    input: bam = f"{ALIGNED_DIR}/aligned.sorted.bam"
    output: f"{ALIGNED_DIR}/aligned.sorted.bam.validation.txt"
    shell: 
        "gatk ValidateSamFile -I {input.bam} -MODE SUMMARY > {output}"

rule duplicates:
    input: bam = f"{ALIGNED_DIR}/aligned.sorted.bam"
    output:
        metrics = f"{ALIGNED_DIR}/dup_metrics.txt",
        dedup = f"{ALIGNED_DIR}/dedup.bam"
    shell: 
        "gatk MarkDuplicates -I {input.bam} -O {output.dedup} -M {output.metrics}"

rule indexing:
    input: bam = f"{ALIGNED_DIR}/dedup.bam"
    output:
        f"{ALIGNED_DIR}/dedup.bam.bai"
    shell: 
        "samtools index {input.bam}"

rule variant_calling:
    input: 
        bam = f"{ALIGNED_DIR}/dedup.bam",
        fasta = f"{RAW_DIR}/reference.fasta"
    output:
        vcf = f"{VARIANT_DIR}/raw_variants.vcf"
    shell: 
        "gatk HaplotypeCaller -R {input.fasta} -I {input.bam} -O {output.vcf}"

rule filtering:
    input: 
        fasta = f"{RAW_DIR}/reference.fasta",
        vcf = f"{VARIANT_DIR}/raw_variants.vcf"
    output:
        vcf = f"{VARIANT_DIR}/filtered_variants.vcf"
    shell: 
        """gatk VariantFiltration -R {input.fasta} -V {input.vcf} -O {output.vcf} --filter-expression "QD < 2.0 || FS > 60.0" --filter-name FILTER"""

rule genbank:
    output:
        gbk = f"{SNPEFF_DATA_DIR}/genes.gbk"
    shell: 
        "efetch -db nucleotide -id {REF_ID} -format genbank > {output.gbk}"

rule custom_snpEff:
    input:
        fasta = f"{RAW_DIR}/reference.fasta",
        gbk = f"{SNPEFF_DATA_DIR}/genes.gbk"
    output:
        config = f"{SNPEFF_DIR}/snpEff.config"
    shell:
        """
        cat <<EOF > {output.config}
        reference_db.genome : reference_db
        reference_db.fa : $(readlink -f {input.fasta})
        reference_db.genbank : $(readlink -f {input.gbk})
        EOF
        """

rule snpEff_database:
    input: config = f"{SNPEFF_DIR}/snpEff.config"
    output: f"{SNPEFF_DIR}/snpEff_reference_db.txt"
    shell:
        "snpEff build -c {input.config} -genbank -v -noCheckProtein reference_db > {output}"

