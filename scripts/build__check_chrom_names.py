#!/usr/bin/env python3
import argparse
import subprocess

"""
Checks that files match the reference genome (ie chromosome names are not mismatched).
Usage: ./gene_mismatch_checker.py <reference FAI file> <...other files (GFF,BED,BAM,VCF) >

    How it works:
        1. Read the reference FAI file and store a map of chromosome names.
        2. Read the other files, and for each, report mismatches of chromosome names.
"""


# ===== MAIN =====
def main():
    # Parse cli args
    parser = argparse.ArgumentParser(
        description="Check that all files match the reference (i.e., no chromosome name mismatches)"
    )
    parser.add_argument("ref", help="The reference FAI file")
    parser.add_argument("files", nargs="*", help="Other files (GFF,BED,BAM,VCF)")

    args = parser.parse_args()

    if not args.ref or not args.ref.endswith(".fai"):
        print(f"Invalid reference file '{args.ref}' (must be a .fai file)")
        exit(1)

    # Read reference file chromsomes
    ref_chroms = fai_extract_chromosomes(args.ref)
    # print(f"\nReference chromosomes from FASTA index:\n{', '.join(ref_chroms)}\n")

    # Read all other files, checking for mismatches against the reference
    n_issues = 0
    for file in args.files:
        chroms = None
        if file.endswith(".gff"):
            chroms = gff_extract_chromosomes(file)
        elif file.endswith(".vcf"):
            chroms = vcf_extract_chromosomes(file)
        elif file.endswith(".bed"):
            chroms = bed_extract_chromosomes(file)
        elif file.endswith(".bam"):
            chroms = bam_extract_chromosomes(file)
        else:
            print(f"Unable to parse '{file}': unrecognised format")
            exit(1)

        chrom_issues = check_mismatches(ref_chroms, list(chroms))
        if len(chrom_issues) > 0:
            print(
                f"Chromosome mismatch{'es' if len(chrom_issues) > 1 else ''} in '{basename(file)}': {', '.join(chrom_issues)}"
            )
            n_issues += 1
    if n_issues == 0:
        exit(0)
    else:
        exit(1)


# Works for fasta.fai or fastq.fai files.
# Returns a dictionary of chromosome names:
# {<chromosome_name>: (length, offset),...}
def fai_extract_chromosomes(file: str) -> dict[str, tuple[int, int]]:
    chroms = {}
    with open(file, "r") as f:
        for line in f:
            fields = line.rstrip().split()
            if len(fields) < 1:
                continue
            feature_len = int(fields[1]) if len(fields) > 1 else None
            feature_offset = int(fields[2]) if len(fields) > 2 else None
            chroms[fields[0]] = (feature_len, feature_offset)
    return chroms


# Works for .gff files.
# Returns a dictionary of chromosome names:
# {<chromosome_name>: (lowest_start_idx, highest_end_idx),...}
def gff_extract_chromosomes(file: str) -> dict[str, tuple[int, int]]:
    chroms = {}
    with open(file, "r") as f:
        for line in f:
            sline = line.rstrip()
            if sline.startswith("#"):
                pass
            else:
                fields = sline.split("\t")
                chrom_name = fields[0]
                feature_start = int(fields[3])
                feature_end = int(fields[4])

                # add chrom to map
                if chrom_name not in chroms:
                    chroms[chrom_name] = (feature_start, feature_end)

                # update feature_start if smaller
                if feature_start < chroms[chrom_name][0]:
                    chroms[chrom_name] = (
                        feature_start,
                        chroms[chrom_name][1],
                    )

                # update feature_end if larger
                if feature_end > chroms[chrom_name][1]:
                    chroms[chrom_name] = (
                        chroms[chrom_name][0],
                        feature_end,
                    )
    return chroms


# Works for .gff files.
# Returns a dictionary of chromosome names:
# {<chromosome_name>: (lowest_start_idx, highest_end_idx),...}
def bed_extract_chromosomes(file: str) -> dict[str, tuple[int, int]]:
    chroms = {}
    with open(file, "r") as f:
        for line in f:
            sline = line.rstrip()
            if sline.startswith("#"):
                pass
                # print("comment", sline)
            else:
                fields = sline.split()
                chrom_name = fields[0]
                feature_start = int(fields[1])
                feature_end = int(fields[2])

                # add chrom to map
                if chrom_name not in chroms:
                    chroms[chrom_name] = (feature_start, feature_end)

                # update feature_start if smaller
                if feature_start < chroms[chrom_name][0]:
                    chroms[chrom_name] = (
                        feature_start,
                        chroms[chrom_name][1],
                    )

                # update feature_end if larger
                if feature_end > chroms[chrom_name][1]:
                    chroms[chrom_name] = (
                        chroms[chrom_name][0],
                        feature_end,
                    )
    return chroms


# Works for .vcf files.
# Returns a dictionary of chromosome names:
# {<chromosome_name>: [pos1, pos2, ..., posX],...}
def vcf_extract_chromosomes(file: str) -> dict[str, list[int]]:
    chroms = {}
    with open(file, "r") as f:
        for line in f:
            sline = line.rstrip()
            if sline.startswith("#"):
                pass
            else:
                fields = sline.split("\t")
                chrom_name = fields[0]
                pos = int(fields[1])

                # add chrom to map
                if chrom_name not in chroms:
                    chroms[chrom_name] = []

                chroms[chrom_name].append(pos)
    return chroms


def bam_extract_chromosomes(file: str) -> dict[str, list[int]]:
    chroms = {}

    cmd = [
        "samtools",
        "view",
        "-H",
        file,
    ]

    result = subprocess.run(cmd, capture_output=True, text=True, check=True)
    for line in result.stdout.splitlines():
        if line.startswith("@SQ"):
            fields = line.split("\t")
            chrom_name = fields[1].split(":")[1]
            chrom_len = int(fields[2].split(":")[1])
            chroms[chrom_name] = (0, chrom_len)
    return chroms


def check_mismatches(refs: dict[str, tuple[int, int]], chroms: list[str]):
    chrom_issues: list[str] = []
    for c in chroms:
        if c not in refs:
            chrom_issues.append(c)
    return chrom_issues


# Helpers
def basename(filename: str) -> str:
    f = filename
    while f.find("/") != -1:
        f = f[f.find("/") + 1 :]
    return f


if __name__ == "__main__":
    main()
