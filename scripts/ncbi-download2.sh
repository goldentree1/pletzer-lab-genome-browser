#!/bin/bash

NCBI_ACCESSION_ID="$1"
OUTPUT_DIR="$2"
DATA_DIR="tmp/build/$NCBI_ACCESSION_ID"
NCBI_DATA_DIR="$DATA_DIR/$NCBI_ACCESSION_ID/ncbi_dataset/data/$NCBI_ACCESSION_ID"

# === MAIN ===

# Download + extract dataset into 'tmp'
mkdir -p "$DATA_DIR"
curl -L \
 -H "Accept: application/zip" \
 -o "$DATA_DIR/$NCBI_ACCESSION_ID.zip" \
 "https://api.ncbi.nlm.nih.gov/datasets/v2/genome/accession/$NCBI_ACCESSION_ID/download?include_annotation_type=GENOME_FASTA&include_annotation_type=GENOME_GFF&include_annotation_type=RNA_FASTA&include_annotation_type=CDS_FASTA&include_annotation_type=PROT_FASTA&include_annotation_type=SEQUENCE_REPORT&hydrated=FULLY_HYDRATED" \
 > /dev/null
unzip "$DATA_DIR/$NCBI_ACCESSION_ID.zip" -d "$DATA_DIR/$NCBI_ACCESSION_ID"

REFSEQ_FILE=$(ls "$NCBI_DATA_DIR" | grep "$NCBI_ACCESSION_ID" | grep ".fna")
GENES_FILE=$(ls "$NCBI_DATA_DIR" | grep "genomic.gff")

# We should really do md5sum checking here too... comes with each NCBI download.

mkdir -p "$DATA_DIR"
cp "$NCBI_DATA_DIR/$REFSEQ_FILE" "$DATA_DIR/refseq.fna"
cp "$NCBI_DATA_DIR/$GENES_FILE" "$DATA_DIR/genomic.gff"

mkdir -p "$OUTPUT_DIR/$NCBI_ACCESSION_ID"
mv "$DATA_DIR/refseq.fna" "$OUTPUT_DIR/$NCBI_ACCESSION_ID/refseq.fna"
mv "$DATA_DIR/genomic.gff" "$OUTPUT_DIR/$NCBI_ACCESSION_ID/genes.gff"
