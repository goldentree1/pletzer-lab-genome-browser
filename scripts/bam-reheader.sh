#!/bin/bash

# Usage:
# reheader-bam.sh <input.bam> <chrom_name> <new_chrom_name>
#
# Example:
# ./scripts/reheader-bam.sh "public/data/LESB58/Pa_BF_s1Aligned.sortedByCoord.out.bam" "chromosome" "NC_011770.1"

samtools view -H "$1" > "$1".sam
sed "s/SN:$2/SN:$3/" "$1".sam > "$1"_FIXED.sam
mv "$1" "$1.ORIGINAL"
samtools reheader "$1"_FIXED.sam "$1.ORIGINAL" > "$1"
rm "$1".sam "$1"_FIXED.sam
