#!/bin/bash

# Define base folder
FOLDER_DIR="/mnt/c/Users/jenhy/Documents/research/ELM-RNAseq/Data/24_Days_RNA_seq"

# Define input/output folders
RAW_DIR="$FOLDER_DIR/raw_data"
TRIMMED_DIR="$FOLDER_DIR/trimmed_data"

RAW_QC_DIR="$FOLDER_DIR/fastqc_pretrim"
TRIMMED_QC_DIR="$FOLDER_DIR/fastqc_posttrim"

# Create output directories if they don't exist
mkdir -p "$RAW_QC_DIR"
mkdir -p "$TRIMMED_QC_DIR"

echo "🔬 Running FastQC on raw paired-end FASTQ files..."

# Run FastQC on raw paired-end reads (R1 and R2)
for FILE in "$RAW_DIR"/*_R[12]_001.fastq.gz; do
  echo "➡️  Processing $(basename "$FILE")"
  docker run --rm \
    -v "$RAW_DIR:/data" \
    -v "$RAW_QC_DIR:/output" \
    pegi3s/fastqc "/data/$(basename "$FILE")" -o /output
done

echo "✅ Pre-trim FastQC complete."
echo

echo "🔬 Running FastQC on trimmed paired-end FASTQ files..."

# Run FastQC on trimmed paired-end reads only (skip unpaired)
for FILE in "$TRIMMED_DIR"/trimmed_*_paired.fastq.gz; do
  echo "➡️  Processing $(basename "$FILE")"
  docker run --rm \
    -v "$TRIMMED_DIR:/data" \
    -v "$TRIMMED_QC_DIR:/output" \
    pegi3s/fastqc "/data/$(basename "$FILE")" -o /output
done

echo "✅ Post-trim FastQC complete."
