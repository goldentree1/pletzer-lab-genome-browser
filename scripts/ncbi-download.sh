#!/bin/bash

NCBI_ACCESSION_ID="$1"
OUTPUT_DIR="$2"
if [[ -z "$OUTPUT_DIR" ]]; then
    OUTPUT_DIR="$NCBI_ACCESSION_ID"
fi

TMP_DIR="tmp/build/$NCBI_ACCESSION_ID"
NCBI_DATA_DIR="$TMP_DIR/$NCBI_ACCESSION_ID/ncbi_dataset/data/$NCBI_ACCESSION_ID"

# === MAIN ===
#
if [[ -d "$TMP_DIR" ]]; then
  rm -rf "$TMP_DIR"
fi

# Download + extract dataset into 'tmp'
mkdir -p "$TMP_DIR"
curl -L \
 -H "Accept: application/zip" \
 -o "$TMP_DIR/$NCBI_ACCESSION_ID.zip" \
 "https://api.ncbi.nlm.nih.gov/datasets/v2/genome/accession/$NCBI_ACCESSION_ID/download?include_annotation_type=GENOME_FASTA&include_annotation_type=GENOME_GFF&include_annotation_type=RNA_FASTA&include_annotation_type=CDS_FASTA&include_annotation_type=PROT_FASTA&include_annotation_type=SEQUENCE_REPORT&hydrated=FULLY_HYDRATED" \
 > /dev/null
unzip "$TMP_DIR/$NCBI_ACCESSION_ID.zip" -d "$TMP_DIR/$NCBI_ACCESSION_ID"

REFSEQ_FILE=$(ls "$NCBI_DATA_DIR" | grep "$NCBI_ACCESSION_ID" | grep ".fna")
GENES_FILE=$(ls "$NCBI_DATA_DIR" | grep "genomic.gff")

# We should really do md5sum checking here too... comes with each NCBI download.

mkdir -p "$TMP_DIR"
cp "$NCBI_DATA_DIR/$REFSEQ_FILE" "$TMP_DIR/refseq.fna"
cp "$NCBI_DATA_DIR/$GENES_FILE" "$TMP_DIR/genomic.gff"

mkdir -p "$OUTPUT_DIR"
mv "$TMP_DIR/refseq.fna" "$OUTPUT_DIR/refseq.fna"
mv "$TMP_DIR/genomic.gff" "$OUTPUT_DIR/genes.gff"
