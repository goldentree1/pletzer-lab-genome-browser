#!/bin/bash

# Usage:
# bam-to-bw.sh <input.bam> <output.bw>
#
# Example:
# ./scripts/bam-to-bw.sh "public/data/LESB58/Pa_BF_s1Aligned.sortedByCoord.out.bam.reheader.bam" "bigwig.bw"

samtools index "$1"
bamCoverage --binSize 5 -b "$1" -o "$2"
