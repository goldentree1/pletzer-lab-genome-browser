# Pletzer Lab Genome Browser
Genome browser website and related scripting utilities.
Made for Pletzer Lab by Elliott Brown.

## **Setup**
On first use, you must setup the dependencies for this project.
To check what you're missing, run:

```bash
scripts/deps-check.sh
```

Then install the listed missing packages.


Next, create & activate the [conda](https://docs.conda.io/projects/conda/en/latest/user-guide/install/index.html) environment:
```bash
conda env create --name plgb --file requirements.conda 
conda activate plgb
```

## **User Guide**

To rebuild the website with new data, follow the steps below:

## **Setup**
On first use, you must setup the dependencies for this project. Check what you're missing:
```bash
scripts/deps-check.sh
```
Then install the listed missing packages.

Next create the [conda](https://docs.conda.io/projects/conda/en/latest/user-guide/install/index.html) environment:
```bash
conda env create --name plgb --file requirements.yaml 
conda activate plgb
```
---

1. **Change into the root directory of this project & activate Conda:**

   ```bash
   cd pletzer-lab-genome-browser/
   ```
---

2. **Create a directory structure like so:**

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
          ├── ...
    ```

---

3. **Run the build script on your data directory:**

    ```bash
    ./scripts/build /path/to/data/ --n-threads=10 --bin-size=10 --yes
    ```

    If the data is malformed, the script will abort and list errors.

    #### Common Errors:
    - Chromosome names are inconsistent:
        A common issue will be that the names of chromosomes are not consistent across files. For example:
        ```bash
        ./scripts/build /path/to/data/
        ```
        ```text
        Checking for errors...
        [OK] GCF_000006765.1 (P.aeruginosa PA01)
        [FAIL] GCF_000013465.1 (S.aureus USA300LAC)
            - Mismatch in 'BF_SA_BF.1.bam': [ CP000255.1 ] not in reference: [ NC_007793.1, NC_007790.1 ]
            - Mismatch in 'BF_SA_BF.2.bam': [ CP000255.1 ] not in reference: [ NC_007793.1, NC_007790.1 ]
        Aborting due to errors.
        ```
    
        If we know that "CP000255.1" chromosome maps to "NC_007793.1", we can fix this by renaming the chromosome in the BAM file to match the reference sequence using `./scripts/bam-reheader`:
        ```bash
        ./scripts/bam-reheader /path/to/BF_SA_BF.1.bam "CP000255.1" "NC_007793.1"
        # After this, there will be two files:
        #  - BF_SA_BF.1.bam.ORIGINAL (the original file with mismatched chromosomes)
        #  - BF_SA_BF.1.bam (the fixed file)
        ```
    ---

4. To preview your website, run the following command:
    ```bash
    npm start
    ```

5. Once happy, copy the entire "dist/" folder to your web server or hosting provider.
   For Wordpress, go to your [Wordpress Admin Dashboard](https://pletzerlab.com/wp-admin), and upload the "dist/" folder using 'WP File Manager' in the sidebar. Once uploaded, rename "dist" to 'pletzer-lab-genome-browser', and it will be available at [https://pletzerlab.com/pletzer-lab-genome-browser](https://pletzerlab.com/pletzer-lab-genome-browser).

## Developer Guide

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
