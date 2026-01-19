#!/bin/bash

rm -rf ./test/tmp
cp -r ./test/test_files ./test/tmp

# get stdout
cmd_stdout=$(scripts/build.sh ./test/tmp | tee /dev/tty)

# Strip ANSI colors just in case
clean_stdout=$(echo "$cmd_stdout" | sed 's/\x1B\[[0-9;]*[mK]//g')

# expected output
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

# Check each line exists in output
failed=0
while IFS= read -r line; do
    # skip empty lines
    [[ -z "$line" ]] && continue
    # check if line exists
    if ! grep -Fxq "$line" <<<"$clean_stdout"; then
        echo "❌ Missing line: $line"
        failed=1
    fi
done <<<"$expected"

if [ "$failed" -eq 0 ]; then
    echo -e "\e[32mPASSED!\e[0m (scripts/build.sh ./test/tmp)"
else
    echo -e "\e[31mFAILED\e[0m (scripts/build.sh ./test/tmp)"
fi
