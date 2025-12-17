#!/bin/bash

# USAGE:
# extract_names_from_bam.sh < src/static/DEMO_DATA/INfixed.bam

samtools view -H | grep @SQ | cut -f2,3
