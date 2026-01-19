#!/bin/bash

main(){
    # First arg is DATA_DIR. Otherwise default to ./data
    DATA_DIR="${1:-./data}"
    BIN_SIZE=10

    # err check: selected dir must exist
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
    # - a "reads" directory (reads) with subfolders for each condition
    # - reads/<condition>/* must contain at least one .bam file
    # - all chromosome names in .bam and .gff files must match reference sequence
    should_exit=false
    echo
    echo "Checking for errors..."
    for genome_dir in "$DATA_DIR"/*; do
        errors=()
        refseq_file="$genome_dir/refseq.fasta"
        genes_file="$genome_dir/genes.gff"
        reads_dir="$genome_dir/reads"

        # err check: refseq.fasta exists (exit: no further checks possible if no refseq)
        if [[ ! -s "$refseq_file" ]]; then
        errors+=("'refseq.fasta' empty or missing — cannot check further")
        printf "\033[0;31m[FAIL]\033[0m \033[0;34m%s\033[0m\n" "$(basename "$genome_dir")"
        for e in "${errors[@]}"; do
            printf "    - %s\n" "$e"
        done
        continue
        fi

        # err check: genes.gff exists
        [[ -f "$genes_file"  ]] || errors+=("'genes.gff' missing")

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

    echo -e "\033[0;32mAll OK!\033[0m"
    echo

    # ================================================
    #     ===   2: Prepare data for JBrowse    ===
    # ================================================
    # If the script gets here, all is well with the data and structure.
    # We can now prepare the data for JBrowse!
    echo "Preparing data for Pletzer Lab Genome Browser..."
    for genome_dir in "$DATA_DIR"/*; do
        refseq_file="$genome_dir/refseq.fasta"
        genes_file="$genome_dir/genes.gff"
        reads_dir="$genome_dir/reads"
        printf "Preparing \033[0;34m%s\033[0m:\n" "$(basename "$genome_dir")"
        echo "Preparing genes file '$(basename "$genes_file")'..."
        jbrowse_prepare_gff "$genes_file"

        # Check and convert BAMs into BigWigs for JBrowse
        for condition_dir in "$reads_dir"/*; do
            [[ -d "$condition_dir" ]] || continue

        bam_files=()
        bw_files=()
        cpm_bw_files=()

        for bam_file in "$condition_dir"/*.bam; do
            [[ -f "$bam_file" ]] || continue
            bam_files+=("$bam_file")

            bam_name=$(basename "$bam_file")
            bw_file="$condition_dir/${bam_name%.bam}.bw"
            cpm_bw_file="$condition_dir/${bam_name%.bam}.cpm.bw"

            echo "Processing BAM file '$bam_name' (this will take a few minutes)..."
            samtools index "$bam_file"
            bamCoverage -b "$bam_file" -o "$bw_file" --binSize "$BIN_SIZE"
            bamCoverage -b "$bam_file" -o "$cpm_bw_file" --normalizeUsing CPM --binSize "$BIN_SIZE"

            bw_files+=("$bw_file")
            cpm_bw_files+=("$cpm_bw_file")
        done

        n_files=${#bam_files[@]}
        if (( n_files == 0 )); then
            echo "No BAM files found in '$condition_dir', skipping."
            continue
        elif (( n_files == 1 )); then
            # Only one BAM: just copy individual BigWigs to “average” names
            avg_bw="$condition_dir/$(basename "$condition_dir").average.bw"
            cpm_avg_bw="$condition_dir/$(basename "$condition_dir").average.cpm.bw"
            echo "Only one BAM file found in '$condition_dir', which will be treated as the average."
            cp "${bw_files[0]}" "$avg_bw"
            cp "${cpm_bw_files[0]}" "$cpm_avg_bw"
        else
            # Multiple BAMs -> merge BAMs and create merged BigWigs
            merged_bam="$condition_dir/$(basename "$condition_dir").merged.bam"
            avg_bw="$condition_dir/$(basename "$condition_dir").average.bw"
            cpm_avg_bw="$condition_dir/$(basename "$condition_dir").average.cpm.bw"

            echo "Merging BAM files into '$merged_bam'"
            samtools merge -f "$merged_bam" "${bam_files[@]}"
            samtools index "$merged_bam"

            echo "Generating merged BigWig"
            bamCoverage -b "$merged_bam" -o "$avg_bw" --binSize "$BIN_SIZE"
            echo "Generating CPM-normalized merged BigWig"
            bamCoverage -b "$merged_bam" -o "$cpm_avg_bw" --normalizeUsing CPM --binSize "$BIN_SIZE"
        fi

        done

    done

}

# Gzip, index and make chrom.sizes a FASTA reference sequence for JBrowse
jbrowse_prepare_fasta(){
    FASTA_FILE=$1
    bgzip -fk "$FASTA_FILE"
    samtools faidx "$FASTA_FILE.gz"
    cut -f1,2 < "$FASTA_FILE.gz.fai" > "$FASTA_FILE.chrom.sizes"
}

# Converts/fixes GFF format and creates indexes for JBrowse
jbrowse_prepare_gff(){
    GFF_FILE=$1

    # Use AGAT to sort and standardise the GFF/GTF file into standard GFF3 format
    agat config --expose --tabix > /dev/null # Make sure there's a config
    agat_convert_sp_gxf2gxf.pl \
      --gff "$GFF_FILE" \
      -o "$GFF_FILE.agat" > /dev/null 2>&1

    # sort chromosomes for tabix
    {
      grep '^#' "$GFF_FILE.agat"
      LC_ALL=C sort -t $'\t' -k1,1 -k4,4n <(grep -v '^#' "$GFF_FILE.agat")
    } > "$GFF_FILE.sorted"

    # remove 'region' types (these just span the entire genome - not useful)
    awk -F'\t' '$3 != "region"' < "$GFF_FILE.sorted" > "$GFF_FILE.sorted.noregion.gff"

    # gzip it and create tabix index
    bgzip -fk "$GFF_FILE.sorted.noregion.gff" -o "$GFF_FILE.sorted.noregion.gff.gz"
    tabix -p gff "$GFF_FILE.sorted.noregion.gff.gz"
}

main "$1"
