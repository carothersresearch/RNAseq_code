#!/bin/bash

# Path to your FASTQ files and adapter file
FOLDER_DIR="/mnt/c/Users/jenhy/Documents/research/ELM-RNAseq/Data/24_Days_RNA_seq"
DATA_DIR="$FOLDER_DIR/raw_data"
ADAPTER="$FOLDER_DIR/TruSeq3-PE.fa"
TRIM_DIR="$FOLDER_DIR/trimmed_data"


# Create output folder if it doesn't exist
mkdir -p "$TRIM_DIR"

# Loop over all R1 files
for R1 in "$DATA_DIR"/*_R1_001.fastq.gz; do
  # Get matching R2
  R2="${R1/_R1_/_R2_}"

  # Extract sample name prefix
  BASE=$(basename "$R1" _R1_001.fastq.gz)

  echo "🔧 Processing $BASE..."

  # Run Trimmomatic using Docker
  docker run --rm \
    -v "$DATA_DIR:/data" \
    -v "$TRIM_DIR:/trim" \
    -v "$FOLDER_DIR:/adapters" \
    staphb/trimmomatic trimmomatic PE -phred33 \
    "/data/${BASE}_R1_001.fastq.gz" \
    "/data/${BASE}_R2_001.fastq.gz" \
    "/trim/trimmed_${BASE}_R1_paired.fastq.gz" \
    "/trim/trimmed_${BASE}_R1_unpaired.fastq.gz" \
    "/trim/trimmed_${BASE}_R2_paired.fastq.gz" \
    "/trim/trimmed_${BASE}_R2_unpaired.fastq.gz" \
    ILLUMINACLIP:/adapters/TruSeq3-PE.fa:2:30:10 \
    LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36
done

echo "✅ Trimming complete. Output saved in: $TRIM_DIR"
