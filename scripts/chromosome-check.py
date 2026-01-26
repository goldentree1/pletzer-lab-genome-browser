#!/usr/bin/env python3
import argparse
import subprocess

"""
Checks that files match the reference genome (ie chromosome names are not mismatched).
Usage: ./gene_mismatch_checker.py <reference FASTA/FAI file> <...other files (GFF,BED,BAM,VCF) >

    How it works:
        1. Read the reference FAI file and store a map of chromosome names.
        2. Read the other files, and for each, report mismatches of chromosome names.
"""


# ===== MAIN =====
def main():
    try:
        # Parse cli args
        parser = argparse.ArgumentParser(
            description="Check that all files match the reference FASTA file (first FASTA/FAI file)"
        )
        parser.add_argument(
            "files", nargs="*", help="Files (FA,FNA,FASTA,FAI,GFF,BED,BAM,VCF)"
        )
        parser.add_argument(
            "-v", "--verbose", action="store_true", help="Enable verbose output"
        )

        args = parser.parse_args()

        ref_chroms_file = None
        ref_chroms = None
        for file in args.files:
            if (
                file.endswith(".fai")
                or file.endswith(".fasta")
                or file.endswith(".fa")
                or file.endswith(".fna")
            ):
                ref_chroms_file = file
                ref_chroms = (
                    fai_extract_chromosomes(file)
                    if file.endswith(".fai")
                    else fasta_extract_chromosomes(file)
                )
                break

        if args.verbose and ref_chroms is not None and ref_chroms_file is not None:
            print(f"--- (refseq) '{basename(ref_chroms_file)}': {ref_chroms}")

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
            elif (
                file.endswith(".fasta")
                or file.endswith(".fa")
                or file.endswith(".fna")
                or file.endswith(".fai")
            ):
                continue
            else:
                print(f"Unable to parse '{file}': unrecognised format")
                exit(1)

            if args.verbose:
                print(f"--- '{basename(file)}': {chroms}")

            if ref_chroms is not None and chroms is not None:
                chrom_issues = check_mismatches(ref_chroms, list(chroms))
                if len(chrom_issues) > 0:
                    print(
                        f"Mismatch{'es' if len(chrom_issues) > 1 else ''} in '{basename(file)}': [ {', '.join(chrom_issues)} ] not in reference: [ {', '.join(ref_chroms)} ]"
                    )
                    n_issues += 1
        if args.verbose and ref_chroms and n_issues == 0:
            print("Success! No mismatches found.")

        exit(0 if n_issues == 0 else 1)

    except Exception as e:
        print(str(e))
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


def fasta_extract_chromosomes(file: str) -> dict[str, tuple[int, int]]:
    chroms = {}

    cmd = ["samtools", "faidx", file, "-o", "-"]

    result = subprocess.run(cmd, capture_output=True, text=True, check=True)
    for line in result.stdout.splitlines():
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
            if sline.startswith("#") or not sline:
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
