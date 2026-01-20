#!/bin/bash
#

ALL_TESTS_PASSED=true

print_passed_test(){
    ARG_TEST_DESCRIPTION="$1"
    ARG_FAILED="$2"
    if [ "$ARG_FAILED" == false ]; then
        echo -e "\n\e[32m--------------------\e[0m"
        echo -e "\n\e[32m   PASSED TEST!\e[0m"
        echo -e "  Passed: $1"
        echo -e "\n\e[32m--------------------\e[0m"
    else
        echo -e "\e[31mTEST FAILED\e[0m (scripts/plgb-build.sh ./tmp)"
    fi
}

# Compare cleaned stdout against expected output line by line
# Usage: check_output "$clean_stdout" "$expected"
# Returns: 0 if all lines found, 1 if any line missing
check_output() {
    local actual="$1"
    local expected="$2"
    local failed=false

    while IFS= read -r line; do
        [[ -z "$line" ]] && continue
        if ! grep -Fxq "$line" <<<"$actual"; then
            echo "❌ Missing line: $line"
            failed=true
        fi
    done <<<"$expected"

    # Return 0 if success, 1 if failed
    if $failed; then
        ALL_TESTS_PASSED=false
        return 1
    else
        return 0
    fi
}

# Clean up previous test files
rm -rf ./tmp
mkdir -p ./tmp/

# ===============================================================
# TEST #1 -- invalid directory provided, which should be rejected
# by the program with error message.

# run cmd and capture output
cmd_stdout=$(scripts/plgb-build.sh -y ./NOT/A/DIR | tee /dev/tty)
clean_stdout=$(echo "$cmd_stdout" | sed 's/\x1B\[[0-9;]*[mK]//g')

# expected output:
read -r -d '' expected <<'EOF'
Directory './NOT/A/DIR' does not exist.
Usage: scripts/plgb-build.sh [options] <data_directory>
Options:
  -y|--yes       Automatically answer 'yes' to all prompts
  -b|--bin-size  Bin size for BigWig files (default: 10)
  -h|--help      Display this help message
EOF

# Check & print success/failure message
if check_output "$clean_stdout" "$expected"; then
    print_passed_test "Invalid directory throws error." false
else
    print_passed_test "Invalid directory throws error." true
fi

# ===============================================================
# TEST #2 -- Valid and invalid test files from 'test/test-files/'
# Three folders are tested. 'simple-invalid' is a minimal, valid directory structure
# with required files, which should pass. 'simple-invalid-empty' is an empty directory
# which should throw an error about no 'refseq.fasta' file. The last being a directory
# with multiple different errors that should be thrown. See expected output below.

# setup test files in a tmp dir
cp -r ./examples/invalid1/ ./tmp/invalid1

# run cmd and capture output
cmd_stdout=$(scripts/plgb-build.sh -y ./tmp/invalid1 | tee /dev/tty)
clean_stdout=$(echo "$cmd_stdout" | sed 's/\x1B\[[0-9;]*[mK]//g')

# expected output:
read -r -d '' expected <<'EOF'
[FAIL] simple-invalid
    - Chromosome mismatch in 'genes.gff': chr5
    - Condition 'empty-sample' has no BAM files
    - Chromosome mismatch in 'sample1.1.bam': notchrom1
    - Chromosome mismatch in 'sample1.2.bam': badchrom2
    - Chromosome mismatch in 'sample1.3.bam': badchr3
    - BAM file 'sample1.badnum.bam' does not have a valid .<number>.bam extension
    - BAM file 'sample2.unnumbered.bam' does not have a valid .<number>.bam extension
[FAIL] simple-invalid-empty
    - 'refseq.fasta' empty or missing — cannot check further
[OK] simple-valid
EOF

# Check & print success/failure message
if check_output "$clean_stdout" "$expected"; then
    print_passed_test "Valid and invalid file structures are handled appropriately." false
else
    print_passed_test "Valid and invalid file structures are handled appropriately." true
fi

# setup test files in a tmp dir
cp -r ./examples/valid1 ./tmp/valid1

# run cmd and capture output

cmd_stdout=$(scripts/plgb-build.sh -y ./tmp/valid1 | tee /dev/tty)
clean_stdout=$(echo "$cmd_stdout" | sed 's/\x1B\[[0-9;]*[mK]//g')

# ===============================================================
# TEST #3 -- Valid test files from 'test/test-files/valid'. This should
# actually run the entire script since there should not be errors.

# TODO make test 3:
# - check that a certain list of files all exist
# - could do some hashsum and compare
for file in ./tmp/simple-valid/reads/sample1/*.bw; do
    echo "Checking file '$(basename "$file")'"
    if [[ -f "$file" ]]; then
        bigWigToBedGraph "$file" "$(dirname "$file")/$(basename "$file" .bw).bedGraph"
    fi
done

if $ALL_TESTS_PASSED; then
    echo -e "\n\e[32m==============================================\e[0m"
    echo -e "\n\n\e[32m    All 'scripts/plgb-build.sh' tests passed!\e[0m"
    echo -e "\n\n\e[32m==============================================\e[0m"
else
    echo -e "\n\e[31m==============================================\e[0m"
    echo -e "\n\n\e[31m    TESTS FAILED FOR 'scripts/plgb-build.sh'!\e[0m"
    echo -e "\n\n\e[31m==============================================\e[0m"
fi
