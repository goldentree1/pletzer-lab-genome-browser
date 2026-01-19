#!/bin/bash

# Usage:
# reheader-bam.sh <input.bam> <chrom_name> <new_chrom_name>
#
# Example:
# ./scripts/reheader-bam.sh "public/data/LESB58/Pa_BF_s1Aligned.sortedByCoord.out.bam" "chromosome" "NC_011770.1"

samtools view -H "$1" > TMPSamHeader.sam
sed "s/SN:$2/SN:$3/" TMPSamHeader.sam > TMPSamHeader_FIXED.sam
mv "$1" "$1.ORIGINAL"
samtools reheader TMPSamHeader_FIXED.sam "$1.ORIGINAL" > "$1"
rm TMPSamHeader.sam TMPSamHeader_FIXED.sam
