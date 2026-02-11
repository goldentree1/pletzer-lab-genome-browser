#!/bin/bash

# Usage:
# bam-to-bw.sh <input.bam> <output.bw>
#
# Example:
# ./scripts/bam-to-bw.sh "public/data/LESB58/Pa_BF_s1Aligned.sortedByCoord.out.bam.reheader.bam" "bigwig.bw"

samtools index "$1"
bamCoverage --binSize 10 -b "$1" -o "$2"
bamCoverage \
  -b "$1" \
  -o "$2.cpm.bw" \
  --normalizeUsing CPM \
  --binSize 10
# ^ sam liked these sizes.
