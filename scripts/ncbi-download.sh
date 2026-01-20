#!/bin/bash

NCBI_ACCESSION_ID="$1"
DATA_DIR="public/data/$NCBI_ACCESSION_ID"
NCBI_DATA_DIR="$DATA_DIR/$NCBI_ACCESSION_ID/ncbi_dataset/data/$NCBI_ACCESSION_ID"

main(){
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

    jbrowse_prepare_fasta "$DATA_DIR/refseq.fna"
    jbrowse_prepare_gff "$DATA_DIR/genomic.gff"

    cd "$DATA_DIR"

    # create config for JBrowse
    jbrowse add-assembly refseq.fna.gz --name "$NCBI_ACCESSION_ID" --load inPlace --type bgzipFasta
    jbrowse add-track genomic.gff.sorted.noregion.gff.gz -a "$NCBI_ACCESSION_ID" --load inPlace -n "Genes" --trackId GFF3GeneTrack
    jbrowse text-index --attributes Name,old_locus_tag,locus_tag --exclude CDS,exon --force
}

jbrowse_prepare_fasta(){
    echo "Preparing FASTA genome reference sequence..."

    FASTA_FILE=$1
    echo "Compressing..."
    bgzip -fk "$FASTA_FILE"

    echo "Creating indexes..."
    samtools faidx "$FASTA_FILE.gz"

    echo "Extracting chromosome sizes..."
    cut -f1,2 < "$FASTA_FILE.gz.fai" > "$FASTA_FILE.chrom.sizes"

    echo "FASTA preparation complete!"
}

jbrowse_prepare_gff(){
    # Converts/fixes GFF format and creates indexes.
    # This allows it to be visualised in JBrowse

    GFF_FILE=$1

    remove_regions() {
        awk -F'\t' '$3 != "region"'
    }

    # GFF/GTF files can have strange file formats.
    # AGAT is a tool built to standardise and sort GFF/GTF files into standard GFF3 format.
    echo "Normalising + sorting GTF/GFF file using AGAT..."
    agat config --expose --tabix
    # agat_convert_sp_gxf2gxf.pl --gff "$GFF_FILE" -o "$GFF_FILE.sorted" > /dev/null

    agat_convert_sp_gxf2gxf.pl \
      --gff "$GFF_FILE" \
      -o "$GFF_FILE.agat" > /dev/null

    # Sort chromosomes for tabix cmd
    {
      grep '^#' "$GFF_FILE.agat"
      LC_ALL=C sort -t $'\t' -k1,1 -k4,4n <(grep -v '^#' "$GFF_FILE.agat")
    } > "$GFF_FILE.sorted"

    # Remove 'region' type genes (these usually just span the entire genome for bacterium - not helpful)
    echo "Removing regions..."
    remove_regions < "$GFF_FILE.sorted" > "$GFF_FILE.sorted.noregion.gff"

    # TODO genome re-naming... or perhaps just JBrowse aliases will work?

    # Jbrowse requires compression with bgzip
    echo "Compressing with bgzip..."
    bgzip -fk "$GFF_FILE.sorted.noregion.gff" -o "$GFF_FILE.sorted.noregion.gff.gz"

    # Jbrowse requires tabix indexes file
    echo "Creating indexes..."
    tabix -p gff "$GFF_FILE.sorted.noregion.gff.gz"

    echo "GFF preparation complete!"

}
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

jbrowse_prepare_fasta "$DATA_DIR/refseq.fna"
jbrowse_prepare_gff "$DATA_DIR/genomic.gff"

cd "$DATA_DIR"

# create config for JBrowse
jbrowse add-assembly refseq.fna.gz --name "$NCBI_ACCESSION_ID" --load inPlace --type bgzipFasta
jbrowse add-track genomic.gff.sorted.noregion.gff.gz -a "$NCBI_ACCESSION_ID" --load inPlace -n "Genes" --trackId GFF3GeneTrack
jbrowse text-index --attributes Name,old_locus_tag,locus_tag --exclude CDS,exon --force
