#!/bin/bash

# Usage: scripts/build.sh [options] <DIRECTORY>
# Example: scripts/build.sh -y --bin-size 50 "path/to/dir"

main() {
    # user-controllable variables via CLI flags + args
    BIN_SIZE=10
    PROMPT_USER_TO_CONTINUE=true
    DATA_DIR="./data"

    # parse flags + args
    while [[ $# -gt 0 ]]; do
        case "$1" in
            -y|--yes)
                PROMPT_USER_TO_CONTINUE=false
                shift
                ;;
            -b|--bin-size)
                BIN_SIZE="$2"
                shift 2
                ;;
            --bin-size=*)
                BIN_SIZE="${1#*=}"
                shift
                ;;
            -h|--help)
                print_help_message
                exit 0
                ;;
            --) # end of flags
                shift
                break
                ;;
            -*)
                echo "Unknown option: $1"
                usage
                exit 1
                ;;
            *)
                DATA_DIR="$1"
                shift
                ;;
        esac
    done

    # err check: arg for dir must exist
    if [[ ! -d "$DATA_DIR" ]]; then
        echo -e "\033[0;31mDirectory '$DATA_DIR' does not exist.\033[0m"
        print_help_message
        exit 1
    fi

    # List the genomes to be processed
    echo
    echo "Found the following genomes:"
    for genome_dir in "$DATA_DIR"/*; do
        echo -e "  - \033[0;34m$(basename "$genome_dir")\033[0m"
    done

    # Prompt user to continue
    if [[ "$PROMPT_USER_TO_CONTINUE" == true ]]; then
        read -r -p "Continue? (Y/n) " reply
        if [[ ! -z "$reply" && ! "$reply" =~ ^[Yy]$ ]]; then
            echo -e "\033[0;31mAborted.\033[0m"
            exit 0
        fi
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
    sleep 0.3
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
        fi

        # err check: genes.gff exists
        [[ -f "$genes_file"  ]] || errors+=("'genes.gff' missing")

        # err check: chromosome names match reference sequence in genes.gff
        if [[ -f "$refseq_file"  ]] && [[ -f "$genes_file" ]]; then
            jbrowse_prepare_fasta "$refseq_file" # need to prep indexes first to check names
            chrom_name_mismatches=$(scripts/build__check_chrom_names.py "$refseq_file.gz.fai" "$genes_file" 2>&1)
            rc=$?
            if (( rc != 0 )); then
                while IFS= read -r line; do
                    errors+=("$line")
                done <<< "$chrom_name_mismatches"
            fi
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

                # 1. At least one BAM file exists
                if [[ ${#bam_files[@]} -eq 0 || ! -f "${bam_files[0]}" ]]; then
                    errors+=("Condition '$condition_name' has no BAM files")
                    continue
                fi

                 # 2. Check that every BAM file has .<number>.bam extension
                for bam in "${bam_files[@]}"; do
                    bam_name=$(basename "$bam")

                    if [[ "$bam_name" != *.bam ]]; then
                        continue
                    fi

                    if [[ ! "$bam_name" =~ \.[0-9]+\.bam$ ]]; then
                        errors+=("BAM file '$bam_name' does not have a valid .<number>.bam extension")
                        continue
                    fi

                    # 3. Check that BAM chromosome names match refseq.fasta
                    # (e.g., refseq with 'chr' and BAM with 'chrom1' is invalid)
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
            rm "$refseq_file".* 2>/dev/null
        fi
    done

    # abort if errors
    if $should_exit; then
        printf "\033[0;31mAborting due to errors.\n\033[0m"
        exit 1
    fi

    echo -e "\033[0;32mAll OK!\n\033[0m"
    sleep 0.5

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
        printf "\033[0;34m%s\033[0m:\n" "$(basename "$genome_dir")"
        echo "Preparing genes file '$(basename "$genes_file")'..."
        jbrowse_prepare_gff "$genes_file"

        # Check and convert BAMs into BigWigs for JBrowse
        for condition_dir in "$reads_dir"/*; do
            [[ -d "$condition_dir" ]] || continue

            bam_files=()
            bw_files=()
            cpm_bw_files=()

            echo "Processing BAM files (reads data). It can take a few minutes for each one to complete..."
            for bam_file in "$condition_dir"/*.bam; do
                [[ -f "$bam_file" ]] || continue
                bam_files+=("$bam_file")

                bam_name=$(basename "$bam_file")
                bw_file="$condition_dir/${bam_name%.bam}.bw"
                cpm_bw_file="$condition_dir/${bam_name%.bam}.cpm.bw"

                echo "Processing BAM file '$bam_name' into BigWig..."
                samtools index "$bam_file"
                bamCoverage -b "$bam_file" -o "$bw_file" --binSize "$BIN_SIZE" > /dev/null 2>&1
                bamCoverage -b "$bam_file" -o "$cpm_bw_file" --normalizeUsing CPM --binSize "$BIN_SIZE" > /dev/null 2>&1

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
                # Multiple BAMs -> merge BAMs and create merged BigWigs + averages
                merged_bam="$condition_dir/$(basename "$condition_dir").merged.bam"
                avg_bw="$condition_dir/$(basename "$condition_dir").average.bw"
                cpm_avg_bw="$condition_dir/$(basename "$condition_dir").average.cpm.bw"

                echo "Detected multiple replicates. Merging BAM files into '$(basename $merged_bam)...'"
                samtools merge -f "$merged_bam" "${bam_files[@]}"
                samtools index "$merged_bam"

                echo "Generating merged BigWig"
                bamCoverage -b "$merged_bam" -o "$avg_bw" --binSize "$BIN_SIZE" > /dev/null 2>&1
                echo "Generating CPM-normalized merged BigWig"
                bamCoverage -b "$merged_bam" -o "$cpm_avg_bw" --normalizeUsing CPM --binSize "$BIN_SIZE" > /dev/null 2>&1
            fi
        done
    done
}

print_help_message() {
    echo "Usage: $0 [options] <data_directory>"
    echo "Options:"
    echo "  -y|--yes       Automatically answer 'yes' to all prompts"
    echo "  -b|--bin-size  Bin size for BigWig files (default: 10)"
    echo "  -h|--help      Display this help message"
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

main "$@"
