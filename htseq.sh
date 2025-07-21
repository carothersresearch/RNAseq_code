#!/bin/bash

# ===== Configurable paths =====
BASE_DIR="/mnt/c/Users/jenhy/Documents/research/ELM-RNAseq/Data/24_Days_RNA_seq"
TRIMMED_DIR="$BASE_DIR/trimmed_data"
GFF_FILE="/data/genome_reference/ELM_ecoli_genome.gff"
OUTPUT_DIR="$BASE_DIR/counts"

# ===== Create output directory if missing =====
mkdir -p "$OUTPUT_DIR"

# ===== Loop through each SAM file =====
for SAM_FILE in "$TRIMMED_DIR"/*.sam; do
    BASENAME=$(basename "$SAM_FILE" .sam)
    OUTPUT_FILE="$OUTPUT_DIR/counts_${BASENAME}.txt"

    if [ ! -f "$OUTPUT_FILE" ]; then
        echo "Processing: ${BASENAME}.sam"
        docker run --rm -v "$BASE_DIR:/data" pegi3s/htseq bash -c \
        "htseq-count -t CDS -i locus_tag /data/trimmed_data/${BASENAME}.sam $GFF_FILE > /data/counts/counts_${BASENAME}.txt"
    else
        echo "Skipped: ${BASENAME}.sam (counts file already exists)"
    fi
done


