#!/bin/bash

# Usage: scripts/build.sh [options] <DIRECTORY>
# Example: scripts/build.sh -y --bin-size 50 "path/to/dir"

main() {
    # user-controllable variables via CLI flags + args
    DATA_DIR="${1%/}"
    BIN_SIZE=10
    PROMPT_TO_CONTINUE=true
    REBUILD_DATA=true
    REBUILD_WEBSITE=true

    # parse flags + args
    while [[ $# -gt 0 ]]; do
        case "$1" in
            -y|--yes)
                PROMPT_TO_CONTINUE=false
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
        genome_dir="${genome_dir%/}"  # remove trailing slash from each genome_dir
        [ -d "$genome_dir" ] || continue
        echo -e "  - \033[0;34m$(basename "$genome_dir")\033[0m"
    done

    # Prompt user to continue
    if [[ "$PROMPT_TO_CONTINUE" == true ]]; then
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
    for genome_dir in "$DATA_DIR"/*; do
        genome_dir="${genome_dir%/}"  # remove trailing slash from each genome_dir
        [ -d "$genome_dir" ] || continue

        errors=()
        refseq_file="$genome_dir/refseq.fasta"
        genes_file="$genome_dir/genes.gff"
        reads_dir="$genome_dir/reads"

        # err check: refseq.fasta exists
        if [[ ! -s "$refseq_file" ]]; then
            errors+=("'refseq.fasta' empty or missing.")
        fi

        # err check: genes.gff exists
        [[ -s "$genes_file"  ]] || errors+=("'genes.gff' empty or missing.")

        # err check: chromosome names match reference sequence in genes.gff
        if [[ -f "$refseq_file"  ]] && [[ -f "$genes_file" ]]; then
            jbrowse_prepare_fasta "$refseq_file" # need to prep indexes first to check names
            chrom_name_mismatches=$(scripts/chromosomes-match-checker.py "$refseq_file.gz.fai" "$genes_file" 2>&1)
            rc=$?
            if (( rc != 0 )); then
                while IFS= read -r line; do
                    errors+=("$line")
                done <<< "$chrom_name_mismatches"
            fi
        fi

        # err check: reads directory exists
        if [[ ! -d "$reads_dir" ]]; then
            errors+=("Reads directory '$reads_dir' does not exist.")
        fi

        # err check:
        # For each reads/* directory:
        # 1. At least one BAM file exists
        # 2. Check that every BAM file has .<number>.bam extension
        # 3. Check that BAM chromosome names match refseq.fasta
        if [[ -f "$refseq_file"  ]] && [[ -d "$reads_dir" ]]; then
            for condition_dir in "$reads_dir"/*; do
                [[ -d "$condition_dir" ]] || continue  # skip if not a directory
                condition_name=$(basename "$condition_dir")
                bam_files=("$condition_dir"/*.bam)

                # err check: at least one BAM file exists in condition directory
                if [[ ${#bam_files[@]} -eq 0 || ! -f "${bam_files[0]}" ]]; then
                    errors+=("Condition '$condition_name' has no BAM files.")
                    continue
                fi

                 # Check BAM files
                for bam_file in "${bam_files[@]}"; do
                    bam_name=$(basename "$bam_file")

                    # skip merged BAMs - we wont process them
                    if [[ "$bam_name" == *.merged.bam ]] ; then
                        continue
                    fi

                    # err check: must follow .<number>.bam extension
                    if [[ ! "$bam_name" =~ \.[0-9]+\.bam$ ]]; then
                        errors+=("BAM file '$bam_name' does not have a valid .<number>.bam extension")
                        continue
                    fi

                    # err check: must be a file (not dir/symlink)
                    if [[ ! -f "$bam_file" ]]; then
                        errors+=("'$bam_file' is not a file (please check the file path).")
                        continue
                    fi

                    # err check: BAM chromosome names must match refseq.fasta
                    # (e.g., refseq with 'chr' and BAM with 'chrom1' is invalid)
                    bam_mismatches=$(scripts/chromosomes-match-checker.py "$refseq_file.gz.fai" "$bam_file" 2>&1)
                    rc=$?
                    if (( rc != 0 )); then
                        while IFS= read -r line; do
                            errors+=("$line")
                        done <<< "$bam_mismatches"
                    fi
                done
            done
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

    # ================================================
    #     ===   2: Prepare data for JBrowse    ===
    # ================================================
    # If the script gets here, all is well with the data and structure.
    # We can now prepare the data for JBrowse!
    echo "Preparing data for Pletzer Lab Genome Browser..."
    generated_config_files=()
    for genome_dir in "$DATA_DIR"/*; do
        genome_dir="${genome_dir%/}"  # remove trailing slash from each genome_dir
        [ -d "$genome_dir" ] || continue
        refseq_file="$genome_dir/refseq.fasta"
        genes_file="$genome_dir/genes.gff"
        reads_dir="$genome_dir/reads"
        printf "\033[0;34m%s\033[0m:\n" "$(basename "$genome_dir")"
        echo "Preparing genes file '$(basename "$genes_file")'..."
        jbrowse_prepare_gff "$genes_file"

        prev_dir=$(pwd)
        cd "$genome_dir"
        # Now genes is prepared, we can use 'jbrowse text-index' for JBrowse search...
        # which requires making an entire config.json file bc JBrowse sucks.
        echo "Preparing JBrowse text searching indexes..."
        [[ -f config.json ]] && rm config.json # rm old conf
        jbrowse add-assembly "refseq.fasta.gz" --name "asm" --load inPlace --type bgzipFasta
        jbrowse add-track "genes.gff.sorted.noregion.gff.gz" -a "asm" --load inPlace -n "Genes" --trackId GFF3GeneTrack
        jbrowse text-index --attributes Name,old_locus_tag,locus_tag --exclude CDS,exon --force
        cd "$prev_dir"

        # Check and convert BAMs into BigWigs for JBrowse
        # Check and convert BAMs into BigWigs for JBrowse
        coverage_json="["  # start JSON array

        echo "Preparing BAM files... (this may take a few minutes for each file)"
        for condition_dir in "$reads_dir"/*; do
            [[ -d "$condition_dir" ]] || continue

            bam_files=()
            bw_files=()
            cpm_bw_files=()

            for bam_file in "$condition_dir"/*.bam; do
                bam_files+=( "$bam_file" )
            done

            # generate individual BigWigs for each BAM file
            for bam_file in "${bam_files[@]}"; do
                bam_name=$(basename "$bam_file" .bam)
                bw_file="$condition_dir/${bam_name}.bw"
                cpm_bw_file="$condition_dir/${bam_name}.cpm.bw"

                # Generate BigWig and CPM-normalized BigWig
                echo "Processing BAM '$bam_name' into BigWig..."
                samtools index "$bam_file"
                bamCoverage -b "$bam_file" -o "$bw_file" --binSize "$BIN_SIZE" > /dev/null 2>&1
                bamCoverage -b "$bam_file" -o "$cpm_bw_file" --normalizeUsing CPM --binSize "$BIN_SIZE" > /dev/null 2>&1

                # add files to coverage array
                bw_files+=( "\"reads/$(basename "$condition_dir")/$bam_name.bw\"" )
                cpm_bw_files+=( "\"reads/$(basename "$condition_dir")/$bam_name.cpm.bw\"" )
            done

            # Generate averaged BigWigs (if >= 2 replicates)
            avg_bw="$condition_dir/$(basename "$condition_dir").average.bw"
            cpm_avg_bw="$condition_dir/$(basename "$condition_dir").average.cpm.bw"
            if (( n_files >= 2 )); then
                merged_bam="$condition_dir/$(basename "$condition_dir").merged.bam"

                # Calc average and CPM-normalized average from merged BAM
                echo "Detected multiple replicates: creating averaged BigWigs from all replicates..."
                samtools merge -f "$merged_bam" "${bam_files[@]}"
                samtools index "$merged_bam"
                bamCoverage -b "$merged_bam" -o "$avg_bw" --binSize "$BIN_SIZE" > /dev/null 2>&1
                bamCoverage -b "$merged_bam" -o "$cpm_avg_bw" --normalizeUsing CPM --binSize "$BIN_SIZE" > /dev/null 2>&1

                # add files to coverage array
                bw_files=( "\"reads/$(basename "$condition_dir")/$(basename "$condition_dir").average.bw\"" "${bw_files[@]}" )
                cpm_bw_files=( "\"reads/$(basename "$condition_dir")/$(basename "$condition_dir").average.cpm.bw\"" "${cpm_bw_files[@]}" )
            fi


            # add condition samples to array
            coverage_json+=$'\n  ['"$(IFS=,; echo "${bw_files[*]}")"','"$(IFS=,; echo "${cpm_bw_files[*]}")"'],'
        done

        # close coverage JSON
        coverage_json="${coverage_json%,}"
        coverage_json+=$'\n]'

        generated_config_file="$genome_dir/generated-config.json"
        generated_config_files+=("$generated_config_file")

        node <<-EOF
        const fs = require('fs');

        const config = {
            "$(basename "$genome_dir")": {
                    ncbiName: "$(basename "$genome_dir")",
                    dataDir: "/data/$(basename "$genome_dir")",
                    firstRegion: "chr1",
                    trixName: "asm",
                    data: {
                        refSeq: "refseq.fasta.gz",
                        genomic: "genes.gff.sorted.noregion.gff.gz",
                        coverage: $coverage_json
                    },
                    extras: []
                }
            };

            fs.writeFileSync("$generated_config_file", JSON.stringify(config, null, 2));
            console.log("Written config to $generated_config_file");
EOF

    done

    merged_config_file="$DATA_DIR/config.json"

    # Convert bash array to quoted JS array literal
    js_array="[ $(printf '%s\n' "${generated_config_files[@]}" | sed 's/.*/"&"/' | paste -sd, -) ]"

    node <<-EOF
    const fs = require('fs');

    const merged = {};
    const files = $js_array;

    for (const f of files) {
        const data = JSON.parse(fs.readFileSync(f, 'utf-8'));
        Object.assign(merged, data);
    }

    fs.writeFileSync("$merged_config_file", JSON.stringify(merged, null, 2));
    console.log("Written combined config to $merged_config_file");
EOF


    echo "Successfully processed your data!"

    # Prompt user to overwrite current website with new data.
    if [[ "$PROMPT_TO_CONTINUE" == true ]]; then
        prompt_user_to_continue "Re-build Pletzer Lab Genome Browser website now?"
    fi

    replace_website_public_data "$DATA_DIR"
}

prompt_user_to_continue(){
    local message="$1"
    read -r -p "$message (Y/n)" reply
    if [[ ! -z "$reply" && ! "$reply" =~ ^[Yy]$ ]]; then
        echo -e "\033[0;31mAborted.\033[0m"
        exit 0
    fi
}

replace_website_public_data(){
    local dir="$1"
    echo "Rebuilding website from path: '$dir'..."
    cp "$dir/config.json" "src/config.json"
    rm -rf ./public/data/
    mkdir -p ./public/data/
    cp -r "$dir/" ./public/
    echo "Success!"
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
    local FASTA_FILE=$1
    bgzip -fk "$FASTA_FILE"
    samtools faidx "$FASTA_FILE.gz"
    cut -f1,2 < "$FASTA_FILE.gz.fai" > "$FASTA_FILE.chrom.sizes"
}

# Converts/fixes GFF format and creates indexes for JBrowse
jbrowse_prepare_gff(){
    local GFF_FILE=$1

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
