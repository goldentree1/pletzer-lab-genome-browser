# Pletzer JBrowse +
Genome browser website and related scripting utilities.
Made for Pletzer Lab by Elliott Brown.

## Overview

This project uses [JBrowse2 Linear Genome View](https://github.com/GMOD/jbrowse-react-linear-genome-view-vite-demo) to build the genome browser website. Scripts are provided to fetch, transform and prepare the data.

## Usage

### Prepare your files

With a dir structure like so:
GCF_000014625/
├── genes.gff
├── reads/
│   ├── PA14_I/
│   │   ├── PA14_I.1.bam
│   │   ├── PA14_I.2.bam
│   │   └── PA14_I.3.bam
│   └── PA14_Un/
│       ├── PA14_Un.1.bam
│       ├── PA14_Un.2.bam
│       └── PA14_Un.3.bam
└── refseq.fna
```bash

```

### To prepare a genome counts comparison
```bash

# === Make dir and get reference genome data from NCBI ===

mkdir -p public/data/GCF_000006765.1

# Download the reference genome from ncbi
scripts/download_and_prepare_from_ncbi.sh GCF_000006765.1

# === Prepare the BAM files (coverage) ===
cp /path/to/bamfile.bam public/data/GCF_000006765.1/ # copy BAM to directory

# 1. Check that the header names match the .fai file
scripts/extract_bam_names.sh < public/data/GCF_000006765.1/bamfile.bam
# 2. Reheader if needed
scripts/reheader-bam.sh public/data/GCF_000006765.1/bamfile.bam "BAD_HEADER" "NEW_HEADER"
# 3. Prepare BigWig from BAM (for compatibility w/ website)
scripts/bam-to-bw.sh public/data/GCF_000006765.1/bamfile.bam bwfile.bw

# === Build the config ===
# Edit src/config.
```

### Install npm dependencies:
```bash
yarn
```

### Run dev server:
```bash
yarn dev
```

### Build for production:
```bash
yarn build # prod-ready static website, output to ./dist/
```

### Run production build:
```bash
yarn start
```

### Project structure

- #### [scripts](./scripts) directory
Scripts and utilities for processing bioinformatics data.

- #### [src](./src) directory
The raw source code for the website.

- #### [public](./public) directory
Static data (e.g., icons, pre-processed bioinformatics data)

- #### [dist](./dist) directory
After generation of appropriate data and building the website, this directory will contain assets and the webpage running JBrowse. Run this directory as a server (for example, with "npx serve -S ./dist") to make the website available on localhost.






# vite with @jbrowse/react-linear-genome-view

This is a demo of using the linear genome view with vite (see
https://vitejs.dev/)

Vite is a build system that is very fast and becoming more popular, using
esbuild and rollup instead of webpack

This particular demo includes several polyfills that are needed for JBrowse
including the Buffer polyfill

## Demo of `@jbrowse/react-linear-genome-view` with vite

See this app running at https://jbrowse.org/demos/lgv-vite/.

## Usage

Run `yarn` and then `yarn dev` to start a development instance

Run `yarn build` which produces a `build` directory that can be deployed to a
static web server
