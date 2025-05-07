#!/usr/bin/env bash

################################################################################
## Configuration

# The name of the directory in raw_data where the RNAseq files are stored
RNA_DIR="RNAseq"

# The name of the directory in raw_data where the RNAseq files are stored
EM_DIR="EMseq"

# The path to the kallisto index to use for quantification
KALLISTO_INDEX="/media/nunez/SSD-Data/genomes/GRCh38_p14_kallisto.idx"

# The path to the EM-seq genome (prepared with Bismark)
EM_SEQ_GENOME="/media/nunez/SSD-Data/genomes/EMseq_hg38"

################################################################################
## Main script

# Set input/output directory variables

INPUT_DIR="raw_data/${RNA_DIR}"
OUTPUT_DIR="output/${RNA_DIR}"

EM_INPUT_DIR="raw_data/${EM_DIR}"
EM_OUTPUT_DIR="output/${EM_DIR}"

# Run FastQC to perform sequencing quality control checks

# mkdir -p $OUTPUT_DIR/fastqc
# 
# fastqc \
#     -t 48 \
#     -o $OUTPUT_DIR/fastqc \
#     $INPUT_DIR/*.fastq.gz

# Aggregate FastQC reports into a single page using MultiQC

# multiqc \
#     --filename $OUTPUT_DIR/fastqc/multiqc.html \
#     $OUTPUT_DIR/fastqc

# Trim adapters (sequences obtained from NEB #E7760 manual)

# mkdir -p $OUTPUT_DIR/trimmed
# 
# # Loop through the forward reads
# for forward in $INPUT_DIR/*_R1_001.fastq.gz; do
#     # Get the sample name
#     name=$(basename ${forward} _R1_001.fastq.gz)
#     echo "Running cutadapt on ${name}..."
#     # Run cutadapt (trimming software)
#     cutadapt \
#         --cores=0 \
#         -m 1 \
#         --poly-a \
#         -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCA \
#         -A AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT \
#         -o "$OUTPUT_DIR/trimmed/${name}_R1.fastq.gz" \
#         -p "$OUTPUT_DIR/trimmed/${name}_R2.fastq.gz" \
#         "$INPUT_DIR/${name}_R1_001.fastq.gz" \
#         "$INPUT_DIR/${name}_R2_001.fastq.gz"
# done
# 
# # Run FastQC again to make sure trimming worked well
# 
# mkdir -p $OUTPUT_DIR/fastqc-trimmed
# 
# fastqc \
#     -t 48 \
#     -o $OUTPUT_DIR/fastqc-trimmed \
#     $OUTPUT_DIR/trimmed/*.fastq.gz
# 
# multiqc \
#     --filename $OUTPUT_DIR/fastqc-trimmed/multiqc.html \
#     $OUTPUT_DIR/fastqc-trimmed
# 
# # Quantify with kallisto
# 
# for forward in $OUTPUT_DIR/trimmed/*_R1.fastq.gz; do
#     name=$(basename ${forward} _R1.fastq.gz)
#     mkdir -p $OUTPUT_DIR/quant/$name
#     kallisto quant \
#         -t 60 \
#         -i $KALLISTO_INDEX \
#         -o $OUTPUT_DIR/quant/$name \
#         $OUTPUT_DIR/trimmed/${name}_R1.fastq.gz \
#         $OUTPUT_DIR/trimmed/${name}_R2.fastq.gz
# done

# Combine all reads into a single TSV

# mkdir -p $OUTPUT_DIR/combined-reads
# 
# uv run combine_reads.py \
#     $OUTPUT_DIR/quant/ \
#     $OUTPUT_DIR/combined-reads/counts.tsv

# Run DESeq2 on samples

# mkdir -p $OUTPUT_DIR/deseq2
# 
# Rscript differential_gene_expression.r \
#     $OUTPUT_DIR/combined-reads/counts.tsv \
#     $OUTPUT_DIR/deseq2/single.csv
# 
# # Generate volcano plot
# 
# mkdir -p $OUTPUT_DIR/volcano-plot
# 
# uv run volcano_plot.py \
#     $OUTPUT_DIR/deseq2 \
#     $OUTPUT_DIR/volcano-plot

# Run FastQC on EMseq

# mkdir -p $EM_OUTPUT_DIR/fastqc

# fastqc \
#     -t 62 \
#     -o $EM_OUTPUT_DIR/fastqc \
#     $EM_INPUT_DIR/*.fastq.gz

# # Run MultiQC on EMseq

# multiqc \
#     --filename $EM_OUTPUT_DIR/fastqc/multiqc.html \
#     $EM_OUTPUT_DIR/fastqc

# # Trim adapters (sequences obtained from NEB #7120 manual)

# mkdir -p $EM_OUTPUT_DIR/trimmed

# # Loop through the forward reads
# for forward in $EM_INPUT_DIR/*_R1_001.fastq.gz; do
#     # Get the sample name
#     name=$(basename ${forward} _R1_001.fastq.gz)
#     echo "Running cutadapt on ${name}..."
#     # Run cutadapt (trimming software)
#     cutadapt \
#         --cores=0 \
#         -m 1 \
#         --poly-a \
#         -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCA \
#         -A AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT \
#         -o "$EM_OUTPUT_DIR/trimmed/${name}_R1.fastq.gz" \
#         -p "$EM_OUTPUT_DIR/trimmed/${name}_R2.fastq.gz" \
#         "$EM_INPUT_DIR/${name}_R1_001.fastq.gz" \
#         "$EM_INPUT_DIR/${name}_R2_001.fastq.gz"
# done

# # Run FastQC again to make sure trimming worked well

# mkdir -p $EM_OUTPUT_DIR/fastqc-trimmed

# fastqc \
#     -t 62 \
#     -o $EM_OUTPUT_DIR/fastqc-trimmed \
#     $EM_OUTPUT_DIR/trimmed/*.fastq.gz

# multiqc \
#     --filename $EM_OUTPUT_DIR/fastqc-trimmed/multiqc.html \
#     $EM_OUTPUT_DIR/fastqc-trimmed

# Bismark alignment

# mkdir -p $EM_OUTPUT_DIR/aligned

# for forward in $EM_INPUT_DIR/*_R1_001.fastq.gz; do
#    name=$(basename ${forward} _R1_001.fastq.gz)
#    bismark \
#       --bam \
#       --parallel 15 \
#       --genome $EM_SEQ_GENOME \
#       -o $EM_OUTPUT_DIR/aligned \
#       -1 $EM_OUTPUT_DIR/trimmed/${name}_R1.fastq.gz \
#       -2 $EM_OUTPUT_DIR/trimmed/${name}_R2.fastq.gz
# done

# Bismark methylation extraction

mkdir -p $EM_OUTPUT_DIR/methylation

for forward in $EM_INPUT_DIR/*_R1_001.fastq.gz; do
   name=$(basename ${forward} _R1_001.fastq.gz)
   bismark_methylation_extractor \
       --parallel 10 \
       --gzip \
       -o $EM_OUTPUT_DIR/methylation \
       $EM_OUTPUT_DIR/aligned/${name}_R1_bismark_bt2_pe.bam
done

# %% Convert Bismark files to bedGraph

mkdir -p $EM_OUTPUT_DIR/methylation-bedGraph

for forward in $EM_INPUT_DIR/*_R1_001.fastq.gz; do
    name=$(basename ${forward} _R1_001.fastq.gz)
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

# Convert Bismark bedGraph files to DSS input files

mkdir -p $EM_OUTPUT_DIR/methylation-dss-separate

for forward in $EM_INPUT_DIR/*_R1_001.fastq.gz; do
    name=$(basename ${forward} _R1_001.fastq.gz)
    uv run bismark_to_dss.py \
        $EM_OUTPUT_DIR/methylation-bedGraph/${name}_OT.txt.gz.bismark.cov.gz \
        $EM_OUTPUT_DIR/methylation-dss-separate/${name}_OT.txt

    uv run bismark_to_dss.py \
        $EM_OUTPUT_DIR/methylation-bedGraph/${name}_OB.txt.gz.bismark.cov.gz \
        $EM_OUTPUT_DIR/methylation-dss-separate/${name}_OB.txt
done

# %% Aggregate original top and bottom strands

mkdir -p $EM_OUTPUT_DIR/methylation-dss-combined

for forward in $EM_INPUT_DIR/*_R1_001.fastq.gz; do
    name=$(basename ${forward} _R1_001.fastq.gz)
    uv run \
        combine_dss.py \
        $EM_OUTPUT_DIR/methylation-dss-separate/${name}_OT.txt \
        $EM_OUTPUT_DIR/methylation-dss-separate/${name}_OB.txt \
        $EM_OUTPUT_DIR/methylation-dss-combined/${name}.txt
done

# %% Filter low read counts

mkdir -p $EM_OUTPUT_DIR/methylation-dss-combined-filtered

uv run filter_dss.py \
    3 \
    $EM_OUTPUT_DIR/methylation-dss-combined \
    $EM_OUTPUT_DIR/methylation-dss-combined-filtered

# %% Compute average methylation information

mkdir -p $EM_OUTPUT_DIR/average-methylation-info

uv run average_methylation_info.py \
    $EM_OUTPUT_DIR/methylation-dss-combined-filtered \
    > $EM_OUTPUT_DIR/average-methylation-info/info.tsv