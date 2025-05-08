#!/usr/bin/env bash

################################################################################
# %% Configuration

# The name of the directory in raw_data where the RNAseq files are stored
RNA_DIR="RNAseq"

# The name of the directory in raw_data where the RNAseq files are stored
EM_DIR="EMseq"

# The path to the kallisto index to use for quantification
KALLISTO_INDEX="/media/nunez/SSD-Data/genomes/GRCh38_p14_kallisto.idx"

# The path to the EM-seq genome (prepared with Bismark)
EM_SEQ_GENOME="/media/nunez/SSD-Data/genomes/EMseq_hg38"

# The number of cores to use
CORES=56

# The number to pass to Bismark's --parallel options
BISMARK_PARALLEL=10

################################################################################
# %% RNA-seq pre-processing

# %% Set input/output directory variables

RNA_INPUT_DIR="raw_data/${RNA_DIR}"
RNA_OUTPUT_DIR="output/${RNA_DIR}"
RNA_SAMPLE_NAMES=( $(cat metadata/rna.csv | cut -d ',' -f1 | tail -n +2) )

# %% Run FastQC to perform sequencing quality control checks

mkdir -p $RNA_OUTPUT_DIR/fastqc

echo "Running fastqc..."

fastqc \
    -t $CORES \
    -o $RNA_OUTPUT_DIR/fastqc \
    $RNA_INPUT_DIR/*.fastq.gz

multiqc \
    --filename $RNA_OUTPUT_DIR/fastqc/multiqc.html \
    $RNA_OUTPUT_DIR/fastqc

# %% Trim adapters (sequences obtained from NEB #E7760 manual)

mkdir -p $RNA_OUTPUT_DIR/trimmed

for name in "${RNA_SAMPLE_NAMES[@]}"; do
    echo "Running cutadapt on ${name}..."
    cutadapt \
        --cores=0 \
        -m 1 \
        --poly-a \
        -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCA \
        -A AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT \
        -o "$RNA_OUTPUT_DIR/trimmed/${name}_R1.fastq.gz" \
        -p "$RNA_OUTPUT_DIR/trimmed/${name}_R2.fastq.gz" \
        "$RNA_INPUT_DIR/${name}_R1_001.fastq.gz" \
        "$RNA_INPUT_DIR/${name}_R2_001.fastq.gz"
done

# %% Run FastQC again to make sure trimming worked well

mkdir -p $RNA_OUTPUT_DIR/fastqc-trimmed

echo "Running fastqc (again)..."

fastqc \
    -t $CORES \
    -o $RNA_OUTPUT_DIR/fastqc-trimmed \
    $RNA_OUTPUT_DIR/trimmed/*.fastq.gz

multiqc \
    --filename $RNA_OUTPUT_DIR/fastqc-trimmed/multiqc.html \
    $RNA_OUTPUT_DIR/fastqc-trimmed

# %% Quantify with kallisto

for name in "${RNA_SAMPLE_NAMES[@]}"; do
    echo "Running kallisto on ${name}..."
    mkdir -p $RNA_OUTPUT_DIR/quant/$name
    kallisto quant \
        -t $CORES \
        -i $KALLISTO_INDEX \
        -o $RNA_OUTPUT_DIR/quant/$name \
        $RNA_OUTPUT_DIR/trimmed/${name}_R1.fastq.gz \
        $RNA_OUTPUT_DIR/trimmed/${name}_R2.fastq.gz
done

# %% Aggregate transcript abundances into gene counts with tximport

mkdir -p $RNA_OUTPUT_DIR/aggregated-reads

echo "Running tximport..."
Rscript \
    scripts/aggregate_transcripts.r \
    metadata/rna.csv \
    $RNA_OUTPUT_DIR/quant/ \
    $RNA_OUTPUT_DIR/aggregated-reads/counts.csv

# %% Perform differential gene expression test with DESeq2

mkdir -p $RNA_OUTPUT_DIR/deseq2

echo "Running DESeq2..."

Rscript \
    scripts/differential_gene_expression.r \
    metadata/rna.csv \
    metadata/comparisons.csv \
    $RNA_OUTPUT_DIR/aggregated-reads/counts.csv \
    $RNA_OUTPUT_DIR/deseq2

################################################################################
# %% EM-seq pre-processing

EM_INPUT_DIR="raw_data/${EM_DIR}"
EM_OUTPUT_DIR="output/${EM_DIR}"
EM_SAMPLE_NAMES=( $(cat metadata/em.csv | cut -d ',' -f1 | tail -n +2) )

# %% Run FastQC on EMseq

mkdir -p $EM_OUTPUT_DIR/fastqc

echo "Running fastqc..."

fastqc \
    -t $CORES \
    -o $EM_OUTPUT_DIR/fastqc \
    $EM_INPUT_DIR/*.fastq.gz

multiqc \
    --filename $EM_OUTPUT_DIR/fastqc/multiqc.html \
    $EM_OUTPUT_DIR/fastqc

# %% Trim adapters (sequences obtained from NEB #7120 manual)

mkdir -p $EM_OUTPUT_DIR/trimmed

for name in "${EM_SAMPLE_NAMES[@]}"; do
    echo "Running cutadapt on ${name}..."
    cutadapt \
        --cores=0 \
        -m 1 \
        -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCA \
        -A AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT \
        -o "$EM_OUTPUT_DIR/trimmed/${name}_R1.fastq.gz" \
        -p "$EM_OUTPUT_DIR/trimmed/${name}_R2.fastq.gz" \
        "$EM_INPUT_DIR/${name}_R1_001.fastq.gz" \
        "$EM_INPUT_DIR/${name}_R2_001.fastq.gz"
done

# %% Run FastQC again to make sure trimming worked well

mkdir -p $EM_OUTPUT_DIR/fastqc-trimmed

echo "Running fastqc (again)..."

fastqc \
    -t $CORES \
    -o $EM_OUTPUT_DIR/fastqc-trimmed \
    $EM_OUTPUT_DIR/trimmed/*.fastq.gz

multiqc \
    --filename $EM_OUTPUT_DIR/fastqc-trimmed/multiqc.html \
    $EM_OUTPUT_DIR/fastqc-trimmed

# %% Bismark alignment (very long-running step!)

mkdir -p $EM_OUTPUT_DIR/aligned

for name in "${EM_SAMPLE_NAMES[@]}"; do
   echo "Running bismark on ${name}..."
   bismark \
      --bam \
      --parallel $BISMARK_PARALLEL \
      --genome $EM_SEQ_GENOME \
      -o $EM_OUTPUT_DIR/aligned \
      -1 $EM_OUTPUT_DIR/trimmed/${name}_R1.fastq.gz \
      -2 $EM_OUTPUT_DIR/trimmed/${name}_R2.fastq.gz
done

# %% Bismark methylation extraction

mkdir -p $EM_OUTPUT_DIR/methylation

for name in "${EM_SAMPLE_NAMES[@]}"; do
   echo "Running bismark_methylation_extractor on ${name}..."
   bismark_methylation_extractor \
       --parallel $BISMARK_PARALLEL \
       --gzip \
       -o $EM_OUTPUT_DIR/methylation \
       $EM_OUTPUT_DIR/aligned/${name}_R1_bismark_bt2_pe.bam
done

# %% Convert Bismark files to bedGraph

mkdir -p $EM_OUTPUT_DIR/methylation-bedGraph

for name in "${EM_SAMPLE_NAMES[@]}"; do
    echo "Running bismark2bedGraph on ${name}..."

    # Original top (OT) strand
    bismark2bedGraph \
        --dir $EM_OUTPUT_DIR/methylation-bedGraph/ \
        -o ${name}_OT.txt \
        $EM_OUTPUT_DIR/methylation/CpG_OT_${name}_R1_bismark_bt2_pe.txt.gz

    # Original bottom (OB) strand
    bismark2bedGraph \
        --dir $EM_OUTPUT_DIR/methylation-bedGraph/ \
        -o ${name}_OB.txt \
        $EM_OUTPUT_DIR/methylation/CpG_OB_${name}_R1_bismark_bt2_pe.txt.gz
done

# %% Convert Bismark bedGraph files to DSS input files

mkdir -p $EM_OUTPUT_DIR/methylation-dss-separate

for name in "${EM_SAMPLE_NAMES[@]}"; do
    echo "Running bismark_to_dss.py on ${name}..."

    uv run \
        --project scripts/ \
        scripts/bismark_to_dss.py \
        $EM_OUTPUT_DIR/methylation-bedGraph/${name}_OT.txt.gz.bismark.cov.gz \
        $EM_OUTPUT_DIR/methylation-dss-separate/${name}_OT.tsv

    uv run \
        --project scripts/ \
        scripts/bismark_to_dss.py \
        $EM_OUTPUT_DIR/methylation-bedGraph/${name}_OB.txt.gz.bismark.cov.gz \
        $EM_OUTPUT_DIR/methylation-dss-separate/${name}_OB.tsv
done

# %% Aggregate original top and bottom strands

mkdir -p $EM_OUTPUT_DIR/methylation-dss-combined

for name in "${EM_SAMPLE_NAMES[@]}"; do
    echo "Running combine_dss.py on ${name}..."
    uv run \
        --project scripts/ \
        scripts/combine_dss.py \
        $EM_OUTPUT_DIR/methylation-dss-separate/${name}_OT.tsv \
        $EM_OUTPUT_DIR/methylation-dss-separate/${name}_OB.tsv \
        $EM_OUTPUT_DIR/methylation-dss-combined/${name}.tsv
done

# %% Filter low read counts

mkdir -p $EM_OUTPUT_DIR/methylation-dss-combined-filtered

echo "Running filter_dss.py..."

uv run \
    --project scripts/ \
    scripts/filter_dss.py \
    3 \
    $EM_OUTPUT_DIR/methylation-dss-combined \
    $EM_OUTPUT_DIR/methylation-dss-combined-filtered

# %% Compute average methylation information

mkdir -p $EM_OUTPUT_DIR/average-methylation-info

echo "Running average_methylation_info.py..."

uv run \
    --project scripts/ \
    scripts/average_methylation_info.py \
    $EM_OUTPUT_DIR/methylation-dss-combined-filtered \
    > $EM_OUTPUT_DIR/average-methylation-info/info.tsv

# %% Perform differential methylation test with DSS

mkdir -p $EM_OUTPUT_DIR/dss

echo "Running DSS..."

Rscript \
    scripts/dss.r \
    metadata/rna.csv \
    metadata/comparisons.csv \
    $RNA_OUTPUT_DIR/aggregated-reads/counts.csv \
    $RNA_OUTPUT_DIR/deseq2

################################################################################
## %% Analysis

# %% Generate volcano plot

mkdir -p output/analysis

uv run \
    --project scripts/ \
    scripts/analyze.py \
    metadata/ \
    $RNA_OUTPUT_DIR/aggregated-reads/counts.csv \
    $RNA_OUTPUT_DIR/deseq2 \
    output/analysis
