#!/bin/bash

jbrowse_prepare_fasta(){
    FASTA_FILE=$1
    bgzip -fk "$FASTA_FILE"
    samtools faidx "$FASTA_FILE.gz"
    cut -f1,2 < "$FASTA_FILE.gz.fai" > "$FASTA_FILE.chrom.sizes"
}

# Converts/fixes GFF format and creates indexes for JBrowse
jbrowse_prepare_gff(){

    GFF_FILE=$1

    # GFF/GTF files can have strange file formats.
    # AGAT is a tool built to standardise and sort GFF/GTF files into standard GFF3 format.
    agat config --expose --tabix > /dev/null # Make sure there's a config
    agat_convert_sp_gxf2gxf.pl \
      --gff "$GFF_FILE" \
      -o "$GFF_FILE.agat" > /dev/null

    # Sort chromosomes for tabix cmd
    {
      grep '^#' "$GFF_FILE.agat"
      LC_ALL=C sort -t $'\t' -k1,1 -k4,4n <(grep -v '^#' "$GFF_FILE.agat")
    } > "$GFF_FILE.sorted"

    # Remove 'region' type genes (these usually just span the entire genome for bacterium - not helpful)
    awk -F'\t' '$3 != "region"' < "$GFF_FILE.sorted" > "$GFF_FILE.sorted.noregion.gff"
    bgzip -fk "$GFF_FILE.sorted.noregion.gff" -o "$GFF_FILE.sorted.noregion.gff.gz"
    tabix -p gff "$GFF_FILE.sorted.noregion.gff.gz"
}

# First arg is DATA_DIR. Otherwise default to ./data
DATA_DIR="${1:-./data}"
if [[ ! -d "$DATA_DIR" ]]; then
    echo "Directory '$DATA_DIR' does not exist."
    exit 1
fi

# List genomes to be processed and prompt user to continue
echo
echo "The following genomes will be prepared:"
for genome_dir in "$DATA_DIR"/*; do
    printf "  - \033[0;34m%s\033[0m\n" "$(basename "$genome_dir")"
done
echo
read -r -p "Continue? (Y/n) " reply

# Exit unless reply (Enter) or Y/y
if [[ ! -z "$reply" && ! "$reply" =~ ^[Yy]$ ]]; then
    printf "\033[0;31mAborted.\033[0m\n"
    exit 1
fi

# ========================================
#     ===   1: Check for errors    ===
# ========================================
# Each directory requires:
# - a "refseq.fasta" (reference sequence)
# - a "genes.gff" (gene names)
# - a "reads" directory (reads)
# - directories under reads/* must contain at least one .bam file
# - all chromosome names in .bam and .gff files match reference sequence
should_exit=false
echo
echo "Checking for errors..."
for genome_dir in "$DATA_DIR"/*; do
    errors=()
    refseq_file="$genome_dir/refseq.fasta"
    genes_file="$genome_dir/genes.gff"
    reads_dir="$genome_dir/reads"

    # err check: refseq.fasta exists
    [[ -f "$refseq_file" ]] || errors+=("refseq.fasta missing")

    # err check: genes.gff exists
    [[ -f "$genes_file"  ]] || errors+=("genes.gff missing")

    # err check: chromosome names match reference sequence in genes.gff
    jbrowse_prepare_fasta "$refseq_file" # need to prep indexes first to check names
    chrom_name_mismatches=$(scripts/build__check_chrom_names.py "$refseq_file.gz.fai" "$genes_file" 2>&1)
    rc=$?
    if (( rc != 0 )); then
        while IFS= read -r line; do
            errors+=("$line")
        done <<< "$chrom_name_mismatches"
    fi

    # err check:
    # For each reads/* directory:
    # 1. At least one BAM file exists
    # 2. Check that every BAM file has .<number>.bam extension
    # 3. Check that BAM chromosome names match refseq.fasta
    if [[ -d "$reads_dir" ]]; then
        for condition_dir in "$reads_dir"/*; do
            [[ -d "$condition_dir" ]] || continue  # skip if not a directory
            condition_name=$(basename "$condition_dir")
            bam_files=("$condition_dir"/*.bam)

            # 1. Check at least one BAM file exists
            if [[ ${#bam_files[@]} -eq 0 || ! -f "${bam_files[0]}" ]]; then
                errors+=("Condition '$condition_name' has no BAM files")
                continue
            fi

            # 2. Check each BAM file is OK...
            for bam in "${bam_files[@]}"; do
                bam_name=$(basename "$bam")

                # Skip any file not ending with .bam
                if [[ "$bam_name" != *.bam ]]; then
                    continue
                fi

                # Check each BAM file extension matches .<int>.bam (e.g., "bamfile.1.bam" is valid)
                if [[ ! "$bam_name" =~ \.[0-9]+\.bam$ ]]; then
                    errors+=("BAM file '$bam_name' does not have a valid .<number>.bam extension")
                    continue
                fi

                # 3. Check BAM chromosome names match reference (e.g., refseq with 'chr' and BAM with 'chrom1' is invalid)
                bam_mismatches=$(scripts/build__check_chrom_names.py "$refseq_file.gz.fai" "$bam" 2>&1)
                rc=$?
                if (( rc != 0 )); then
                    while IFS= read -r line; do
                        errors+=("$line")
                    done <<< "$bam_mismatches"
                fi
            done
        done
    else
        errors+=("Reads directory '$reads_dir' is missing")
    fi

    if (( ${#errors[@]} == 0 )); then
        printf "\033[0;32m[OK]\033[0m \033[0;34m%s\033[0m\n" "$(basename "$genome_dir")"
    else
        should_exit=true
        printf "\033[0;31m[FAIL]\033[0m \033[0;34m%s\033[0m\n" "$(basename "$genome_dir")"
        for e in "${errors[@]}"; do
            printf "    - %s\n" "$e"
        done
        # remove generated fasta stuff (keeps OG file)
        rm "$refseq_file".*
    fi
done

# abort if errors
if $should_exit; then
    printf "\033[0;31mAborting due to errors.\033[0m\n"
    exit 1
fi

echo "All OK! Preparing data for Pletzer Lab Genome Browser..."
# ================================================
#     ===   2: Prepare data for JBrowse    ===
# ================================================
# If the script gets here, all is well with the data and structure
# We can now prepare the data for JBrowse!
for genome_dir in "$DATA_DIR"/*; do
    refseq_file="$genome_dir/refseq.fasta"
    genes_file="$genome_dir/genes.gff"
    reads_dir="$genome_dir/reads"

    jbrowse_prepare_gff "$genes_file"
done
