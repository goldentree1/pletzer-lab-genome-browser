#!/bin/bash

NCBI_ACCESSION_ID="$1"

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
    #!/bin/bash

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
      -o "$GFF_FILE.agat"

    {
      grep '^#' "$GFF_FILE.agat"
      LC_ALL=C sort -t $'\t' -k1,1 -k4,4n <(grep -v '^#' "$GFF_FILE.agat")
    } > "$GFF_FILE.sorted"

    # # Sort CHROMOSOMES for tabix
    # grep -v '^#' "$GFF_FILE.agat" \
    # | sort -k1,1 -k4,4n \
    # > "$GFF_FILE.sorted"

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
mkdir -p tmp
curl -L \
 -H "Accept: application/zip" \
 -o "tmp/$NCBI_ACCESSION_ID.zip" \
 "https://api.ncbi.nlm.nih.gov/datasets/v2/genome/accession/$NCBI_ACCESSION_ID/download?include_annotation_type=GENOME_FASTA&include_annotation_type=GENOME_GFF&include_annotation_type=RNA_FASTA&include_annotation_type=CDS_FASTA&include_annotation_type=PROT_FASTA&include_annotation_type=SEQUENCE_REPORT&hydrated=FULLY_HYDRATED" \
 > /dev/null
unzip "tmp/$NCBI_ACCESSION_ID.zip" -d "tmp/$NCBI_ACCESSION_ID"

# We should really do md5sum checking here too... comes with each NCBI download.

DATA_FILES_DIR="tmp/$NCBI_ACCESSION_ID/ncbi_dataset/data/$NCBI_ACCESSION_ID"
REFSEQ_FILE=$(ls "$DATA_FILES_DIR" | grep "$NCBI_ACCESSION_ID" | grep ".fna")
GENES_FILE=$(ls "$DATA_FILES_DIR" | grep "genomic.gff")

mkdir -p "public/data/$NCBI_ACCESSION_ID"
cp "$DATA_FILES_DIR/$REFSEQ_FILE" "public/data/$NCBI_ACCESSION_ID/refseq.fna"
cp "$DATA_FILES_DIR/$GENES_FILE" "public/data/$NCBI_ACCESSION_ID/genomic.gff"

jbrowse_prepare_fasta "public/data/$NCBI_ACCESSION_ID/refseq.fna"
jbrowse_prepare_gff "public/data/$NCBI_ACCESSION_ID/genomic.gff"

cd "public/data/$NCBI_ACCESSION_ID"

# create config for JBrowse
jbrowse add-assembly refseq.fna.gz --name "$NCBI_ACCESSION_ID" --load inPlace --type bgzipFasta
jbrowse add-track genomic.gff.sorted.noregion.gff.gz -a "$NCBI_ACCESSION_ID" --load inPlace -n "Genes" --trackId GFF3GeneTrack
jbrowse text-index
