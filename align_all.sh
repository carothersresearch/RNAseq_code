#!/bin/bash

# Set up paths
FOLDER_DIR="/mnt/c/Users/jenhy/Documents/research/ELM-RNAseq/Data/24_Days_RNA_seq"
TRIM_DIR="$FOLDER_DIR/trimmed_data"
REF_DIR="$FOLDER_DIR/genome_reference"
INDEX_BASE="/ref/ELM_ecoli_genome_index"

# Loop through all R1 paired reads
for R1 in "$TRIM_DIR"/trimmed_*_R1_paired.fastq.gz; do
  # Infer sample name
  BASE=$(basename "$R1" _R1_paired.fastq.gz)
  R2="$TRIM_DIR/${BASE}_R2_paired.fastq.gz"
  SAM_OUT="$TRIM_DIR/${BASE}.sam"

  echo "🔄 Aligning $BASE..."

  docker run --rm \
    -v "$TRIM_DIR:/reads" \
    -v "$REF_DIR:/ref" \
    staphb/bowtie2 \
    bowtie2 -x "$INDEX_BASE" \
    -1 "/reads/$(basename "$R1")" -2 "/reads/$(basename "$R2")" \
    -S "/reads/$(basename "$SAM_OUT")" -p 4

  echo "✅ Done: $BASE.sam created."
done

echo "🎉 All alignments completed."

