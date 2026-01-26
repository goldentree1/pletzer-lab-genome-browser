#!/bin/bash

# Usage: scripts/build.sh [options] <DIRECTORY>
# Example: scripts/build.sh -y --bin-size 50 "path/to/dir"

main() {

    # user-controllable variables
    DATA_DIR="${1%/}"
    BIN_SIZE=10 # Bin size used for bioinformatics calculations with BAM/BigWig files.
    PROMPT_TO_CONTINUE=true # Give user prompt of (Y\n) to continue.
    REBUILD_BIGWIGS=true # always rebuild bigwigs
    SKIP_PROCESSED_BAMS=false # if a .bam already has .bw of same name, skip it.
    SKIP_WEBSITE=false
    CHECK_ONLY=false
    N_THREADS=1

    # parse flags + args
    while [[ $# -gt 0 ]]; do
        case "$1" in
            -y|--yes)
                PROMPT_TO_CONTINUE=false
                shift
                ;;
            -c|--check)
                CHECK_ONLY=true
                shift
                ;;
            --skip-bams)
                REBUILD_BIGWIGS=false
                shift
                ;;
            --skip-processed-bams)
                SKIP_PROCESSED_BAMS=true
                shift
                ;;
            --skip-website)
                SKIP_WEBSITE=true
                shift
                ;;
            --n-threads)
                N_THREADS="$2"
                shift 2
                ;;
            --n-threads=*)
                N_THREADS="${1#*=}"
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
                print_full_help
                exit 0
                ;;
            --) # end of flags
                shift
                break
                ;;
            -*)
                echo "Invalid option: $1"
                print_minimal_help
                exit 1
                ;;
            *)
                DATA_DIR="$1"
                shift
                ;;
        esac
    done

    if [[ "$DATA_DIR" == "" ]]; then
        echo -e "\033[0;31mA directory must be provided.\033[0m"
        print_minimal_help
        exit 1
    fi

    # exit if user provided invalid directory
    if [[ ! -d "$DATA_DIR" ]]  then
        echo -e "\033[0;31mDirectory '$DATA_DIR' does not exist or is not readable.\033[0m"
        print_minimal_help
        exit 1
    fi

    # check for errors + print them
    if ! genome_error_check_routine "$DATA_DIR"; then
        printf "\033[0;31mAborting due to errors.\n\033[0m"
        exit 1
    fi

    echo -e "\033[0;32mAll OK!\033[0m"
    echo

    if [[ "$CHECK_ONLY" == true ]]; then
        exit 0
    fi

    # No errors! Prompt user to continue (unless -y flag)
    if [[ "$PROMPT_TO_CONTINUE" == true ]]; then
        read -r -p "Begin processing data? (Y/n) " reply
        if [[ ! -z "$reply" && ! "$reply" =~ ^[Yy]$ ]]; then
            echo -e "\033[0;31mAborted.\033[0m"
            exit 0
        fi
    fi

    # Data is fine, prepare it for JBrowse!!
    echo "Preparing data for Pletzer Lab Genome Browser..."
    generated_config_files=()
    for genome_dir in "$DATA_DIR"/*; do
        genome_dir="${genome_dir%/}"  # remove trailing slash
        refseq_file="$genome_dir/refseq.fasta"
        genes_file="$genome_dir/genes.gff"
        reads_dir="$genome_dir/reads"

        [ -d "$genome_dir" ] || continue # ignore non-directories
        printf " --- \033[0;34m$(basename "$genome_dir")\033[0m --- \n"

        # Create all necessary files for JBrowse
        genome_data_processing_ritual "$refseq_file" "$genes_file" "$genome_dir" "$N_THREADS" "$BIN_SIZE" "$REBUILD_BIGWIGS" "$SKIP_PROCESSED_BAMS"
    done
    # TODO make this function and catch Javascript errors and throw if so

    echo "Data processing is complete!"

    echo "Replacing site configuration..."
    merged_config_file="$DATA_DIR/config.json"
    merge_json_configs "$merged_config_file" "${generated_config_files[@]}"
    # TODO catch Javascript errors and throw if so

    echo "Overwriting data in 'public/data'..."
    replace_public_data "$DATA_DIR"
    # TODO check for cp errs??
    echo "Removing unnecessary files from 'public/data'..."
    clean_public_data
    echo

    if [[ "$SKIP_WEBSITE" == true ]]; then
        echo "Complete!"
        exit 0
    fi

    # Prompt user to overwrite current website with new data.
    if [[ "$PROMPT_TO_CONTINUE" == true ]]; then
        read -r -p "Re-build Pletzer Lab Genome Browser website now? (Y/n) " reply
        if [[ ! -z "$reply" && ! "$reply" =~ ^[Yy]$ ]]; then
            echo "You can rebuild it manually by running:"
            echo "      npm run build"
            echo -e "\033[0;31mExiting.\033[0m"
            exit 0
        fi
    fi

    # Re-build the website
    echo "Building website to 'dist/'... (this may take a minute or two)"
    npm run build > /dev/null 2>&1
    echo "Complete!"
}


# -------
# HELPERS
# -------
#

print_minimal_help() {
    echo "Usage: $0 [options] <data_directory>"
    echo "Run '$0 --help' for full options."
}

print_full_help() {
    echo "Usage: $0 [options] <data_directory>"
    echo
    echo "Options:"
    echo "  -h, --help                  Display this help message"
    echo "  -y, --yes                   Do not prompt for confirmation"
    echo "  -c, --check                 Check for errors without building"
    echo "      --skip-bams             Skip rebuilding BigWig files from BAMs"
    echo "      --skip-processed-bams   Skip rebuilding BigWig files that already exist"
    echo "      --skip-website          Skip rebuilding the website in dist/"
    echo "      --n-threads <n>         Number of threads to use"
    echo "      --n-threads=<n>         "
    echo "  -b, --bin-size <n>          Bin size for BigWig files (default: 10)"
    echo "      --bin-size=<n>          "
    echo
    echo "Arguments:"
    echo "  <data_directory>              Directory with following structure:"
    echo "                                - refseq.fasta"
    echo "                                - genes.gff"
    echo "                                - reads/<condition>/*.<N>.bam"
    echo "  (Where <condition> is the name of the condition and <N> is the number of the replicate)"
    echo
    echo "Example:"
    echo "  scripts/build.sh --yes --bin-size=50 /path/to/<data_directory>"
}



# Error checker for input files.
# Each directory (genome strain) requires:
# - a "refseq.fasta" (reference sequence)
# - a "genes.gff" (gene names)
# - a "reads" directory (reads) with subfolders for each condition
# - reads/<condition>/* must contain at least one .bam file
# - all chromosome names in .bam and .gff files must match reference sequence
genome_error_check_routine(){
    local data_dir="$1"
    should_exit=false
    echo
    echo "Checking for errors..."
    for genome_dir in "$data_dir"/*; do
        genome_dir="${genome_dir%/}"  # remove trailing slash from each genome_dir
        [ -d "$genome_dir" ] || continue

        errors=()
        local refseq_file="$genome_dir/refseq.fasta"
        local genes_file="$genome_dir/genes.gff"
        local reads_dir="$genome_dir/reads"

        # err check: refseq.fasta exists
        if [[ ! -s "$refseq_file" ]]; then
            errors+=("'refseq.fasta' empty or missing.")
        fi

        # err check: genes.gff exists
        [[ -s "$genes_file"  ]] || errors+=("'genes.gff' empty or missing.")

        # err check: chromosome names match reference sequence in genes.gff
        if [[ -f "$refseq_file"  ]] && [[ -f "$genes_file" ]]; then
            # jbrowse_prepare_fasta "$refseq_file" # need fai indexes to check chromosome names
            chrom_name_mismatches=$(scripts/chromosome-check.py "$refseq_file" "$genes_file" 2>&1)
            ret_code=$?
            if (( ret_code != 0 )); then
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

                condition_name=$(basename "$condition_dir")
                bam_files=("$condition_dir"/*.bam)

                [[ -d "$condition_dir" ]] || continue  # skip non-directories

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

                done

                # err check: BAM chromosome names must match refseq.fasta
                # (e.g., refseq with 'chr' and BAM with 'chrom1' is invalid)
                bam_mismatches=$(scripts/chromosome-check.py "$refseq_file" "${bam_files[@]}")
                rc=$?
                if (( rc != 0 )); then
                    while IFS= read -r line; do
                        errors+=("$line")
                    done <<< "$bam_mismatches"
                fi
            done
        fi

        if (( ${#errors[@]} == 0 )); then
            # no errors! all ok.
            printf "\033[0;32m[OK]\033[0m \033[0;34m%s\033[0m\n" "$(basename "$genome_dir")"
        else
            # errors - show them and exit
            should_exit=true
            rm "$refseq_file".* 2>/dev/null # remove generated fasta stuff (keeps OG file)
            printf "\033[0;31m[FAIL]\033[0m \033[0;34m%s\033[0m\n" "$(basename "$genome_dir")"
            for e in "${errors[@]}"; do
                printf "    - %s\n" "$e"
            done
        fi
    done

    $should_exit && return 1 || return 0
}

genome_data_processing_ritual(){
    local refseq_file="$1"
    local genes_file="$2"
    local genome_dir="$3"
    local n_threads="$4"
    local bin_size="$5"
    local rebuild_bigwigs="$6"
    local skip_processed_bams="$7"

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

    # Prepare reference sequence
    echo "Preparing reference sequence '$(basename "$refseq_file")'..."
    jbrowse_prepare_fasta "$refseq_file"

    # Prepare genes file
    echo "Preparing genes file '$(basename "$genes_file")'..."
    jbrowse_prepare_gff "$genes_file"

    # Now genes is prepared, we can use 'jbrowse text-index' for JBrowse search...
    # which requires making an entire config.json file bc JBrowse sucks.
    echo "Preparing JBrowse text searching indexes..."
    prev_dir=$(pwd)
    cd "$genome_dir"
    [[ -f config.json ]] && rm config.json # rm old conf
    jbrowse add-assembly "$(basename "$refseq_file").gz" --name "asm" --load inPlace --type bgzipFasta
    jbrowse add-track "$(basename "$genes_file").sorted.noregion.gff.gz" -a "asm" --load inPlace -n "Genes" --trackId GFF3GeneTrack
    jbrowse text-index --attributes Name,old_locus_tag,locus_tag --exclude CDS,exon --force
    cd "$prev_dir"


    condition_names_json="["
    # Check and convert BAMs into BigWigs for JBrowse
    coverage_json="["  # start JSON array

    echo "Preparing BAM files... (this can take a few minutes for each file)"
    for condition_dir in "$reads_dir"/*; do

        [[ -d "$condition_dir" ]] || continue # skip non-directories
        condition_names_json+="\"$(basename "$condition_dir")\","

        # collect allll the BAM files
        n_replicates=0
        bam_files=()
        bw_files=()
        cpm_bw_files=()
        # These are for bigwigAverage to use
        raw_bws_for_avg=()
        cpm_bws_for_avg=()

        avg_bw="$condition_dir/$(basename "$condition_dir").average.bw"
        cpm_avg_bw="$condition_dir/$(basename "$condition_dir").average.cpm.bw"
        for bam_file in "$condition_dir"/*.bam; do
            if [[ "$(basename "$bam_file")" == *.merged.bam ]]; then
              continue
            fi
            bam_files+=( "$bam_file" )
            n_replicates=$((n_replicates + 1))
        done

        # generate individual BigWigs for each BAM file
        for bam_file in "${bam_files[@]}"; do
            bam_name=$(basename "$bam_file" .bam)
            bw_file="$condition_dir/${bam_name}.bw"
            cpm_bw_file="$condition_dir/${bam_name}.cpm.bw"

            if [[ "$rebuild_bigwigs" == true ]]; then
                if [[ "$skip_processed_bams" == true ]] && [[ -f "$bw_file" && -f "$cpm_bw_file" ]]; then
                    echo "Skipping BAM '$bam_name' as BigWigs already exist."
                else
                    echo "Processing BAM '$bam_name' into BigWig..."
                    samtools index "$bam_file"
                    bamCoverage --numberOfProcessors "$n_threads" -b "$bam_file" -o "$bw_file" --binSize "$bin_size" > /dev/null 2>&1
                    echo "Processing BAM '$bam_name' into BigWig (with CPM-normalisation)..."
                    bamCoverage --numberOfProcessors "$n_threads" -b "$bam_file" -o "$cpm_bw_file" --normalizeUsing CPM --binSize "$bin_size" > /dev/null 2>&1
                fi
            fi

            # Update JSON strings
            bw_files+=( "\"reads/$(basename "$condition_dir")/$bam_name.bw\"" )
            cpm_bw_files+=( "\"reads/$(basename "$condition_dir")/$bam_name.cpm.bw\"" )

            # Keep clean paths for bigwigAverage
            raw_bws_for_avg+=("$bw_file")
            cpm_bws_for_avg+=("$cpm_bw_file")
        done

        # Generate averaged BigWigs (>= 2 replicates)
        if (( n_replicates >= 2 )); then
            echo "Detected multiple replicates to be averaged..."
            if [[ "$rebuild_bigwigs" == true ]]; then
                if [[ "$skip_processed_bams" == true ]] && [[ -f "$avg_bw" && -f "$cpm_avg_bw" ]]; then
                    echo "Skipping already-processed BigWig averages..."
                else
                    echo "Generating BigWig average..."
                    # Normal average
                    bigwigAverage --numberOfProcessors "$n_threads" --binSize "$bin_size" -b "${raw_bws_for_avg[@]}" -o "$avg_bw" > /dev/null 2>&1
                    # CPMs averaged
                    echo "Generating BigWig average (CPM-normalized)..."
                    bigwigAverage --numberOfProcessors "$n_threads" --binSize "$bin_size" -b "${cpm_bws_for_avg[@]}" -o "$cpm_avg_bw" > /dev/null 2>&1
                fi
            fi
            bw_files=( "\"reads/$(basename "$condition_dir")/$(basename "$condition_dir").average.bw\"" "${bw_files[@]}" )
            cpm_bw_files=( "\"reads/$(basename "$condition_dir")/$(basename "$condition_dir").average.cpm.bw\"" "${cpm_bw_files[@]}" )
        fi

        # disgusting hackery to make a JSON array of BigWigs... in bash
        coverage_json+=$'\n  ['"$(IFS=,; echo "${bw_files[*]}")"','"$(IFS=,; echo "${cpm_bw_files[*]}")"'],'
    done


    coverage_json="${coverage_json%,}"
    coverage_json+=$'\n]'
    condition_names_json+="]"
    generated_config_file="$genome_dir/generated-config.json"
    generated_config_files+=("$generated_config_file")

    first_region=$(head -n 1 "$refseq_file.gz.fai" | awk '{print $1}')
    # Optional sanity check
    if [[ -z "$first_region" ]]; then
        echo "ERROR: Could not determine first region from FAI" >&2
        exit 1
    fi
    echo "Found first display region: $first_region"

    # Evil Javascript hackery to finally build CONFIG!! :O
    # Writes all these bash vars into correct config format for the website's Javascript
    node <<-EOF
    const fs = require('fs');

    const config = {
        "$(basename "$genome_dir")": {
                ncbiName: "$(basename "$genome_dir")",
                firstRegion: "$first_region",
                trixName: "asm",
                data: {
                    refSeq: "refseq.fasta.gz",
                    genomic: "genes.gff.sorted.noregion.gff.gz",
                    coverage: $coverage_json.map((arr) => arr.sort((a, b) => {
                        if (a.endsWith(".cpm.bw") && !b.endsWith(".cpm.bw")) return -1; // a goes first
                        if (!a.endsWith(".cpm.bw") && b.endsWith(".cpm.bw")) return 1;  // b goes first
                        return 0;
                    })),
                    coverage_condition_names: $condition_names_json
                },
                extras: []
            }
        };

        fs.writeFileSync("$generated_config_file", JSON.stringify(config, null, 2));
        console.log("Written config to $generated_config_file");
EOF
}

merge_json_configs() {
    local output_file="$1"
    shift
    local files=("$@")
    local js_array

    # yuck: crafting a JS array in bash
    js_array="[ $(printf '%s\n' "${files[@]}" | sed 's/.*/"&"/' | paste -sd, -) ]"

    # this might be even worse hackery...
    # i gave up and just used javascript to read each file and parse
    node <<-EOF
    const fs = require('fs');

    const merged = {};
    const files = $js_array;

    for (const f of files) {
        const data = JSON.parse(fs.readFileSync(f, 'utf-8'));
        Object.assign(merged, data);
    }

    fs.writeFileSync("$output_file", JSON.stringify(merged, null, 2));
    console.log("Written combined config to $output_file");
EOF
}

replace_public_data(){
    local dir="$1"
    cp "$dir/config.json" "config.json"
    rm -rf ./public/data/
    cp -r "$dir/" ./public/data
}

clean_public_data(){
    local dir="./public/data"
    find "$dir" -type f -name "*.bam" -delete
    find "$dir" -type f -name "*.bam.bai" -delete
    find "$dir" -type f -name "*.bam.ORIGINAL" -delete
    find "$dir" -type f -name "*.gff" -delete
    find "$dir" -type f -name "*.gff.sorted" -delete
    find "$dir" -type f -name "*.agat" -delete
    find "$dir" -type f -name "*.fasta" -delete
    find "$dir" -type f -name "*.fa" -delete
    find "$dir" -type f -name "*.fna" -delete
    find "$dir" -maxdepth 2 -type f -name "*.json" -delete
}

main "$@"
