# Pletzer Lab Genome Browser
Genome browser website and related scripting utilities.
Made for Pletzer Lab by Elliott Brown.

## User Guide

To rebuild the website with new data, follow the steps below:

1. Change into the root directory of this project:

    ```bash
    cd pletzer-lab-genome-browser/
    ```

    <strong>If this is your first time using this project, you will need to setup a few things:</strong>
    
    - First ensure you have the necessary dependencies installed:
      - Node.js
      - npm
      - Python 3
    
    - Create a conda environment and install requirements:
    
      ```bash
      conda create -n plgb python=3.8
      conda activate plgb
      pip install -r requirements.txt
      ```
    
    - Install website dependencies:
    
      ```bash
      npm install
      ```

2. Create a directory structure like so:

    ```text
    data/
    ├── GCF_000014625/
    │   ├── refseq.fna
    │   ├── genes.gff
    │   └── reads/
    │       ├── PA14_I/
    │       │   ├── PA14_I.1.bam
    │       │   ├── PA14_I.2.bam
    │       │   └── PA14_I.3.bam
    │       └── PA14_Un/
    │           ├── PA14_Un.1.bam
    │           ├── PA14_Un.2.bam
    │           └── PA14_Un.3.bam
    └── GCF_000026645.1/
        ├── refseq.fna
        └── ...
    ```

3. Run the build script on your directory:

    ```bash
    scripts/build.sh /path/to/data/
    ```

    If there are errors in the data or your directory structure is incorrect, the script will abort and list the errors it encountered. You will need to fix these errors and try again. 

    A common issue will be that the names of chromosomes are not consistent across files. For example:

    ```bash
    (plgb) user@computer:~/PletzerLabGenomeBrowser$ scripts/build.sh -y /path/to/data/
    ```

## Development

### Install npm dependencies:
```bash
npm install
```

### Run dev server:
```bash
npm run dev
```

### Build for production:
```bash
npm run build # prod-ready static website, output to ./dist/
```

### Run production build:
```bash
npm start
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


### Overview

This project uses [JBrowse2 Linear Genome View](https://github.com/GMOD/jbrowse-react-linear-genome-view-vite-demo) to build the genome browser website. Scripts are provided to fetch, transform and prepare the data.
